#' Estimate facet normals from a regularised facet-centre envelope
#'
#' Reconstructs a continuous triangular envelope from the facet-centre point
#' cloud and estimates one outward unit normal per facet from an area- and
#' distance-weighted average of nearby envelope-face normals.
#'
#' The envelope topology is derived directly from facet positions. Each facet
#' is linked to its six nearest facet centres, only reciprocal links are
#' retained, and triangular three-cliques of that graph are used as candidate
#' envelope faces. Candidate faces are rejected when any edge is longer than
#' 1.4 times the larger local facet-spacing estimate at its endpoints. This
#' prevents long bridging faces across gaps while keeping the construction
#' independent of the optical-neighbour graph used downstream.
#'
#' Before normals are calculated, the envelope vertices undergo five fixed
#' Taubin-style smoothing cycles (lambda = 0.30, mu = -0.31), which suppress
#' high-frequency positional noise while limiting shrinkage of the eye surface.
#'
#' For each facet, local facet spacing is the median distance to its six nearest
#' facet centres (or all available neighbours when fewer than six exist). The
#' `envelope_factor` multiplies this spacing to define the Gaussian distance
#' scale (sigma) used to weight envelope-face normals. Face weights are
#' proportional to face area times `exp(-0.5 * (distance / sigma)^2)` and are
#' truncated at 2.5 sigma. No subsequent neighbour-normal averaging is applied.
#'
#' CV3D currently assumes all coordinates are in micrometres (µm). Consequently
#' the returned support scales are also in µm.
#'
#' @param df Data frame or tibble containing facet centres. Required columns are
#'   `ID`, `x`, `y`, and `z`.
#' @param envelope_factor Positive numeric scalar controlling the spatial scale
#'   of envelope-normal estimation relative to local facet spacing. Values used
#'   by the CV3D UI are 1, 1.25, 1.5, and 2; the default is 1.25.
#' @param verbose Logical; print concise progress messages.
#'
#' @return A tibble with one row per facet and columns `ID`, `norm.x`,
#'   `norm.y`, `norm.z`, `normal_method`, `normal_envelope_factor`,
#'   `normal_support_scale_um`, `normal_weight_cutoff_um`, and
#'   `normal_support_face_count`.
#'
#' @details
#' The triangle-based [get_facet_normals()] estimator remains available as an
#' alternative for method comparison and sensitivity analyses.
#'
#' @seealso [get_facet_normals()]
#'
#' @export
get_facet_normals_envelope <- function(df,
                                       envelope_factor = 1.25,
                                       verbose = FALSE) {
  required <- c("ID", "x", "y", "z")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "get_facet_normals_envelope() is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  envelope_factor <- suppressWarnings(as.numeric(envelope_factor)[1])
  if (!is.finite(envelope_factor) || envelope_factor <= 0) {
    stop("envelope_factor must be a finite number greater than zero.", call. = FALSE)
  }

  work <- tibble::as_tibble(df)
  finite_rows <- is.finite(suppressWarnings(as.numeric(work$x))) &
    is.finite(suppressWarnings(as.numeric(work$y))) &
    is.finite(suppressWarnings(as.numeric(work$z)))
  work <- work[finite_rows, ]
  work$ID <- as.character(work$ID)
  work$x <- as.numeric(work$x)
  work$y <- as.numeric(work$y)
  work$z <- as.numeric(work$z)

  if (nrow(work) < 4) {
    stop("At least four finite facet centres are required.", call. = FALSE)
  }
  if (anyDuplicated(work$ID) > 0) {
    stop("Facet IDs must be unique.", call. = FALSE)
  }

  n <- nrow(work)
  coords <- cbind(work$x, work$y, work$z)
  colnames(coords) <- c("x", "y", "z")

  # Six-nearest-neighbour geometry is used only to reconstruct the envelope and
  # to define its local scale. It is deliberately independent of the 05A
  # optical-neighbour graph.
  k_points <- min(7L, n)
  nearest <- RANN::nn2(data = coords, query = coords, k = k_points)
  nearest_idx <- nearest$nn.idx
  nearest_dist <- nearest$nn.dists
  if (k_points > 1L) {
    nearest_idx <- nearest_idx[, -1, drop = FALSE]
    nearest_dist <- nearest_dist[, -1, drop = FALSE]
  }

  local_spacing <- apply(nearest_dist, 1, stats::median, na.rm = TRUE)
  global_spacing <- stats::median(
    local_spacing[is.finite(local_spacing) & local_spacing > 0],
    na.rm = TRUE
  )
  if (!is.finite(global_spacing) || global_spacing <= 0) {
    stop("Could not estimate local facet spacing from facet-centre positions.", call. = FALSE)
  }
  local_spacing[!is.finite(local_spacing) | local_spacing <= 0] <- global_spacing

  # Mutual six-nearest-neighbour graph.
  reciprocal <- vector("list", n)
  for (i in seq_len(n)) {
    candidates <- unique(as.integer(nearest_idx[i, ]))
    candidates <- candidates[is.finite(candidates) & candidates >= 1L & candidates <= n & candidates != i]
    if (length(candidates) == 0) {
      reciprocal[[i]] <- integer(0)
      next
    }
    keep <- vapply(
      candidates,
      function(j) i %in% as.integer(nearest_idx[j, ]),
      logical(1)
    )
    reciprocal[[i]] <- candidates[keep]
  }

  # Triangular three-cliques of the mutual graph. Construct each once by
  # requiring i < a < b.
  face_list <- vector("list", n)
  face_count <- 0L
  for (i in seq_len(n)) {
    candidates <- sort(reciprocal[[i]][reciprocal[[i]] > i])
    if (length(candidates) < 2) next
    pairs <- utils::combn(candidates, 2)
    for (k in seq_len(ncol(pairs))) {
      a <- pairs[1, k]
      b <- pairs[2, k]
      if (b %in% reciprocal[[a]]) {
        face_count <- face_count + 1L
        face_list[[face_count]] <- c(i, a, b)
      }
    }
  }
  if (face_count < 1L) {
    stop("No triangular faces could be constructed from facet-centre positions.", call. = FALSE)
  }
  faces <- do.call(rbind, face_list[seq_len(face_count)])
  storage.mode(faces) <- "integer"

  # Remove bridging triangles. This reproduces the envelope topology used in
  # CV3D's degradation validation: every face edge must remain within 1.4 times
  # the larger local spacing at its two endpoints.
  edge_ok <- function(i, j) {
    d <- sqrt(rowSums((coords[j, , drop = FALSE] - coords[i, , drop = FALSE])^2))
    limit <- 1.4 * pmax(local_spacing[i], local_spacing[j])
    is.finite(d) & is.finite(limit) & d <= limit
  }
  keep_face <- edge_ok(faces[, 1], faces[, 2]) &
    edge_ok(faces[, 2], faces[, 3]) &
    edge_ok(faces[, 1], faces[, 3])
  faces <- faces[keep_face, , drop = FALSE]
  if (nrow(faces) < 1L) {
    stop("All candidate envelope faces were rejected as bridging faces.", call. = FALSE)
  }

  # Use only edges represented by retained faces for the smoothing adjacency.
  edges <- unique(rbind(
    faces[, c(1, 2), drop = FALSE],
    faces[, c(2, 3), drop = FALSE],
    faces[, c(1, 3), drop = FALSE]
  ))
  storage.mode(edges) <- "integer"
  edge_from <- c(edges[, 1], edges[, 2])
  edge_to <- c(edges[, 2], edges[, 1])
  degree <- tabulate(edge_from, nbins = n)

  neighbour_mean <- function(xyz) {
    out <- xyz
    sums <- rowsum(xyz[edge_to, , drop = FALSE], group = edge_from, reorder = FALSE)
    groups <- suppressWarnings(as.integer(rownames(sums)))
    valid <- is.finite(groups) & groups >= 1L & groups <= n & degree[groups] > 0
    groups <- groups[valid]
    if (length(groups) > 0) {
      out[groups, ] <- sums[valid, , drop = FALSE] / degree[groups]
    }
    out
  }

  # Five low-shrinkage Taubin-style cycles, matching the validation test.
  smoothed <- coords
  for (iter in seq_len(5L)) {
    mean_xyz <- neighbour_mean(smoothed)
    smoothed <- smoothed + 0.30 * (mean_xyz - smoothed)
    mean_xyz <- neighbour_mean(smoothed)
    smoothed <- smoothed - 0.31 * (mean_xyz - smoothed)
  }

  a <- smoothed[faces[, 1], , drop = FALSE]
  b <- smoothed[faces[, 2], , drop = FALSE]
  c <- smoothed[faces[, 3], , drop = FALSE]
  ab <- b - a
  ac <- c - a
  cross <- cbind(
    ab[, 2] * ac[, 3] - ab[, 3] * ac[, 2],
    ab[, 3] * ac[, 1] - ab[, 1] * ac[, 3],
    ab[, 1] * ac[, 2] - ab[, 2] * ac[, 1]
  )
  cross_length <- sqrt(rowSums(cross^2))
  valid_face <- is.finite(cross_length) & cross_length > .Machine$double.eps
  if (!any(valid_face)) {
    stop("All reconstructed envelope faces are degenerate.", call. = FALSE)
  }

  faces <- faces[valid_face, , drop = FALSE]
  a <- a[valid_face, , drop = FALSE]
  b <- b[valid_face, , drop = FALSE]
  c <- c[valid_face, , drop = FALSE]
  cross <- cross[valid_face, , drop = FALSE]
  cross_length <- cross_length[valid_face]

  face_normals <- cross / cross_length
  face_area <- cross_length / 2
  face_centres <- (a + b + c) / 3

  eye_centroid_smoothed <- colMeans(smoothed)
  outward_score <- rowSums(
    face_normals * sweep(face_centres, 2, eye_centroid_smoothed, FUN = "-")
  )
  flip_face <- is.finite(outward_score) & outward_score < 0
  face_normals[flip_face, ] <- -face_normals[flip_face, , drop = FALSE]

  sigma <- local_spacing * envelope_factor
  cutoff <- 2.5 * sigma

  # Query a conservative bounded number of nearby faces rather than allocating
  # the full facet-by-face distance matrix. For the validated support range this
  # is numerically equivalent to the full calculation while scaling to large eyes.
  k_faces <- min(
    nrow(face_centres),
    max(64L, as.integer(ceiling(64 * envelope_factor^2)))
  )
  batch_size <- 5000L
  normal_matrix <- matrix(NA_real_, nrow = n, ncol = 3)
  support_face_count <- integer(n)

  starts <- seq.int(1L, n, by = batch_size)
  for (start in starts) {
    idx <- start:min(n, start + batch_size - 1L)
    face_nn <- RANN::nn2(
      data = face_centres,
      query = coords[idx, , drop = FALSE],
      k = k_faces
    )

    nn_idx <- face_nn$nn.idx
    nn_dist <- face_nn$nn.dists
    row_sigma <- sigma[idx]
    ratio <- nn_dist / row_sigma

    area_matrix <- matrix(
      face_area[as.vector(nn_idx)],
      nrow = nrow(nn_idx), ncol = ncol(nn_idx)
    )
    weights <- area_matrix * exp(-0.5 * ratio^2)
    weights[!is.finite(weights) | !is.finite(ratio) | ratio > 2.5] <- 0

    zero_rows <- which(rowSums(weights) <= 0)
    if (length(zero_rows) > 0) {
      weights[cbind(zero_rows, rep(1L, length(zero_rows)))] <-
        area_matrix[cbind(zero_rows, rep(1L, length(zero_rows)))]
    }

    support_face_count[idx] <- rowSums(weights > 0)

    fx <- matrix(
      face_normals[, 1][as.vector(nn_idx)],
      nrow = nrow(nn_idx), ncol = ncol(nn_idx)
    )
    fy <- matrix(
      face_normals[, 2][as.vector(nn_idx)],
      nrow = nrow(nn_idx), ncol = ncol(nn_idx)
    )
    fz <- matrix(
      face_normals[, 3][as.vector(nn_idx)],
      nrow = nrow(nn_idx), ncol = ncol(nn_idx)
    )

    weighted <- cbind(
      rowSums(weights * fx),
      rowSums(weights * fy),
      rowSums(weights * fz)
    )
    weighted_length <- sqrt(rowSums(weighted^2))
    good <- is.finite(weighted_length) & weighted_length > .Machine$double.eps
    weighted[good, ] <- weighted[good, , drop = FALSE] / weighted_length[good]
    weighted[!good, ] <- NA_real_
    normal_matrix[idx, ] <- weighted
  }

  # Orient the final vectors outwards relative to the original point cloud.
  eye_centroid <- colMeans(coords)
  outward <- rowSums(
    normal_matrix * sweep(coords, 2, eye_centroid, FUN = "-")
  )
  flip <- is.finite(outward) & outward < 0
  normal_matrix[flip, ] <- -normal_matrix[flip, , drop = FALSE]

  if (isTRUE(verbose)) {
    message(
      "Envelope normals: ", n, " facets, ", nrow(face_centres),
      " faces, factor=", format(envelope_factor, trim = TRUE), "."
    )
  }

  tibble::tibble(
    ID = work$ID,
    norm.x = normal_matrix[, 1],
    norm.y = normal_matrix[, 2],
    norm.z = normal_matrix[, 3],
    normal_method = "envelope",
    normal_envelope_factor = envelope_factor,
    normal_support_scale_um = sigma,
    normal_weight_cutoff_um = cutoff,
    normal_support_face_count = support_face_count
  )
}
