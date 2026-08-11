#' Calculate the Angle Between Two 3D Vectors
#'
#' Calculates the smallest angle between two vectors in three-dimensional
#' space from their normalized dot product.
#'
#' @param a Numeric vector of length three.
#' @param b Numeric vector of length three.
#' @param unit Character. Unit of the returned angle. Either `"radians"` or
#'   `"degrees"`. Default: `"radians"`.
#'
#' @return A numeric scalar containing the angle between `a` and `b`. The
#'   returned angle lies between 0 and pi radians, or between 0 and 180
#'   degrees. If either vector has zero length, `NA_real_` is returned.
#'
#' @keywords internal
vector_angle <- function(a, b, unit = c("radians", "degrees")) {
  unit <- match.arg(unit)
  a <- as.numeric(a)
  b <- as.numeric(b)

  if (length(a) != 3 || length(b) != 3 || any(!is.finite(a)) || any(!is.finite(b))) {
    stop("Both vectors must contain three finite numeric values.", call. = FALSE)
  }

  magnitude_a <- sqrt(sum(a^2))
  magnitude_b <- sqrt(sum(b^2))
  if (magnitude_a == 0 || magnitude_b == 0) return(NA_real_)

  cos_angle <- sum(a * b) / (magnitude_a * magnitude_b)
  cos_angle <- max(-1, min(1, cos_angle))
  angle <- acos(cos_angle)

  if (unit == "degrees") angle <- angle * 180 / pi
  angle
}


#' Identify Neighbouring Facets
#'
#' Identifies likely neighbouring facets from their three-dimensional
#' positions. Facet coordinates are projected onto a unit sphere, so
#' neighbourhood relationships are evaluated by angular rather than Cartesian
#' distance.
#'
#' Candidate neighbour links are obtained from a convex-hull triangulation of
#' the normalized facet positions. Links that are long relative to the local
#' angular facet spacing are removed. Facets with too few remaining neighbours
#' are supplemented with their nearest angular neighbours, and the number of
#' neighbours can optionally be limited.
#'
#' @param df A data frame or tibble containing one row per facet.
#' @param x Character. Name of the column containing x coordinates. Default:
#'   `"x"`.
#' @param y Character. Name of the column containing y coordinates. Default:
#'   `"y"`.
#' @param z Character. Name of the column containing z coordinates. Default:
#'   `"z"`.
#' @param id Character. Name of the column containing unique facet identifiers.
#'   Default: `"ID"`.
#' @param center Logical. If `TRUE`, centre the point cloud on its coordinate
#'   centroid before projecting points onto the unit sphere. Default: `TRUE`.
#' @param k_local Integer. Number of nearest angular neighbours used to
#'   estimate the local facet spacing. Default: `6`.
#' @param knn_search Integer. Number of nearest angular neighbours retained as
#'   candidates for supplementing facets with too few triangulation
#'   neighbours. Default: `20`.
#' @param edge_tol Numeric. Relative tolerance for retaining triangulation
#'   edges. An edge is retained when its angular length does not exceed
#'   `(1 + edge_tol)` times the larger local-spacing estimate of its two
#'   endpoints. Default: `0.5`.
#' @param min_neighbours Integer. Minimum number of neighbours targeted by
#'   the nearest-neighbour fallback for sparse facets. Default: `3`.
#' @param max_neighbours Integer or `NULL`. Maximum number of reciprocal
#'   neighbours retained for each facet. If `NULL`, no maximum is applied.
#'   Default: `6`.
#'
#' @return The input data with two additional columns: `neighbours`, containing
#'   the IDs of neighbouring facets separated by `"; "`, and
#'   `number_of_neighbours`, containing the number of neighbours assigned to
#'   each facet. Neighbour relationships are reciprocal: if facet A lists facet
#'   B, facet B also lists facet A.
#'
#' @examples
#' data(cv3d_example_facets)
#' neighbours <- find_neighbours(cv3d_example_facets)
#' head(neighbours[, c("ID", "neighbours", "number_of_neighbours")])
#'
#' @export
find_neighbours <- function(df,
                            x = "x", y = "y", z = "z",
                            id = "ID",
                            center = TRUE,
                            k_local = 6,
                            knn_search = 20,
                            edge_tol = 0.5,
                            min_neighbours = 3,
                            max_neighbours = 6) {
  if (!is.data.frame(df)) stop("df must be a data frame or tibble.", call. = FALSE)
  if (!all(c(x, y, z, id) %in% names(df))) {
    stop("Missing required coordinate or ID columns.", call. = FALSE)
  }
  if (anyDuplicated(df[[id]]) > 0) stop("IDs must be unique.", call. = FALSE)
  if (!is.logical(center) || length(center) != 1 || is.na(center)) {
    stop("center must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(k_local) || length(k_local) != 1 || !is.finite(k_local) ||
      k_local < 1 || k_local != as.integer(k_local)) {
    stop("k_local must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(knn_search) || length(knn_search) != 1 || !is.finite(knn_search) ||
      knn_search < 1 || knn_search != as.integer(knn_search)) {
    stop("knn_search must be a positive integer.", call. = FALSE)
  }
  if (knn_search < k_local) stop("knn_search must be >= k_local.", call. = FALSE)
  if (!is.numeric(edge_tol) || length(edge_tol) != 1 || !is.finite(edge_tol) || edge_tol < 0) {
    stop("edge_tol must be a single finite number >= 0.", call. = FALSE)
  }
  if (!is.numeric(min_neighbours) || length(min_neighbours) != 1 ||
      !is.finite(min_neighbours) || min_neighbours < 0 || min_neighbours != as.integer(min_neighbours)) {
    stop("min_neighbours must be a non-negative integer.", call. = FALSE)
  }
  if (!is.null(max_neighbours) &&
      (!is.numeric(max_neighbours) || length(max_neighbours) != 1 ||
       !is.finite(max_neighbours) || max_neighbours < 1 ||
       max_neighbours != as.integer(max_neighbours))) {
    stop("max_neighbours must be NULL or a positive integer.", call. = FALSE)
  }
  min_neighbours <- as.integer(min_neighbours)
  k_local <- as.integer(k_local)
  knn_search <- as.integer(knn_search)
  if (!is.null(max_neighbours) && min_neighbours > max_neighbours) {
    stop("min_neighbours cannot exceed max_neighbours.", call. = FALSE)
  }

  pts <- as.matrix(df[, c(x, y, z), drop = FALSE])
  storage.mode(pts) <- "double"
  if (nrow(pts) < 4) stop("At least four facet positions are required.", call. = FALSE)
  if (any(!is.finite(pts))) stop("Facet coordinates must be finite.", call. = FALSE)

  n <- nrow(pts)
  if (center) pts <- sweep(pts, 2, colMeans(pts), "-")

  norms <- sqrt(rowSums(pts^2))
  if (any(norms == 0)) stop("Some points have zero norm after centering; cannot normalize.", call. = FALSE)
  u <- pts / norms

  ang_dist <- function(ui, uj) {
    d <- sum(ui * uj)
    d <- max(min(d, 1), -1)
    acos(d)
  }

  knn_eff <- min(as.integer(knn_search), n - 1L)
  k_local_eff <- min(as.integer(k_local), knn_eff)
  nn <- RANN::nn2(u, u, k = knn_eff + 1L)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  dch <- nn$nn.dists[, -1, drop = FALSE]
  dang <- 2 * asin(pmin(dch / 2, 1))

  local_scale <- apply(dang[, seq_len(k_local_eff), drop = FALSE], 1, stats::median, na.rm = TRUE)

  tri <- geometry::convhulln(u, options = "Qt")
  if (is.null(tri) || nrow(tri) == 0) stop("convhulln failed (degenerate point set?).", call. = FALSE)

  edge_pairs <- function(a, b) cbind(pmin(a, b), pmax(a, b))
  e <- rbind(
    edge_pairs(tri[, 1], tri[, 2]),
    edge_pairs(tri[, 2], tri[, 3]),
    edge_pairs(tri[, 3], tri[, 1])
  )
  e <- unique(e)
  colnames(e) <- c("i", "j")

  keep <- logical(nrow(e))
  for (k in seq_len(nrow(e))) {
    i <- e[k, "i"]
    j <- e[k, "j"]
    threshold <- (1 + edge_tol) * max(local_scale[i], local_scale[j])
    keep[k] <- ang_dist(u[i, ], u[j, ]) <= threshold
  }
  e <- e[keep, , drop = FALSE]

  adj <- replicate(n, integer(0), simplify = FALSE)
  if (nrow(e) > 0) {
    for (k in seq_len(nrow(e))) {
      i <- e[k, "i"]
      j <- e[k, "j"]
      adj[[i]] <- c(adj[[i]], j)
      adj[[j]] <- c(adj[[j]], i)
    }
    adj <- lapply(adj, unique)
  }

  minimum_neighbours <- min(min_neighbours, k_local_eff)
  if (!is.null(max_neighbours)) {
    max_neighbours <- as.integer(max_neighbours)
    minimum_neighbours <- min(minimum_neighbours, max_neighbours)
  }

  # Supplement sparse triangulation regions with nearest angular neighbours.
  # Every added relationship is inserted in both directions.
  for (i in seq_len(n)) {
    if (length(adj[[i]]) < minimum_neighbours) {
      cand <- idx[i, ]
      for (j in cand) {
        if (length(adj[[i]]) >= minimum_neighbours) break
        if (j == i || j %in% adj[[i]]) next
        adj[[i]] <- c(adj[[i]], j)
        adj[[j]] <- unique(c(adj[[j]], i))
      }
    }
  }

  if (!is.null(max_neighbours)) {
    edge_length <- function(i, j) ang_dist(u[i, ], u[j, ])

    # Symmetrically prune the longest edges until every facet satisfies the
    # requested maximum. Prefer removals that do not leave either endpoint
    # below the fallback minimum.
    repeat {
      degree <- lengths(adj)
      over <- which(degree > max_neighbours)
      if (length(over) == 0) break

      edge_mat <- do.call(rbind, lapply(over, function(i) {
        if (length(adj[[i]]) == 0) return(NULL)
        cbind(i = i, j = adj[[i]])
      }))
      if (is.null(edge_mat) || nrow(edge_mat) == 0) break
      edge_mat <- cbind(
        i = pmin(edge_mat[, "i"], edge_mat[, "j"]),
        j = pmax(edge_mat[, "i"], edge_mat[, "j"])
      )
      edge_mat <- unique(edge_mat)

      removable <- degree[edge_mat[, "i"]] > minimum_neighbours &
        degree[edge_mat[, "j"]] > minimum_neighbours
      candidates <- edge_mat[removable, , drop = FALSE]
      if (nrow(candidates) == 0) candidates <- edge_mat

      lengths_now <- vapply(seq_len(nrow(candidates)), function(k) {
        edge_length(candidates[k, "i"], candidates[k, "j"])
      }, numeric(1))
      k_remove <- which.max(lengths_now)
      i <- candidates[k_remove, "i"]
      j <- candidates[k_remove, "j"]
      adj[[i]] <- setdiff(adj[[i]], j)
      adj[[j]] <- setdiff(adj[[j]], i)
    }

    # If pruning made a facet sparse, refill from nearest angular neighbours
    # only where both endpoints still have capacity.
    repeat {
      degree <- lengths(adj)
      sparse <- which(degree < minimum_neighbours)
      if (length(sparse) == 0) break
      changed <- FALSE

      for (i in sparse) {
        if (length(adj[[i]]) >= minimum_neighbours ||
            length(adj[[i]]) >= max_neighbours) next
        cand <- idx[i, ]
        for (j in cand) {
          if (j == i || j %in% adj[[i]]) next
          if (length(adj[[j]]) >= max_neighbours) next
          adj[[i]] <- c(adj[[i]], j)
          adj[[j]] <- unique(c(adj[[j]], i))
          changed <- TRUE
          break
        }
      }

      if (!changed) break
    }
  }

  adj <- lapply(adj, function(v) sort(unique(v)))


  ids <- df[[id]]
  df$neighbours <- vapply(adj, function(v) {
    if (length(v) == 0) return("")
    paste(ids[v], collapse = "; ")
  }, character(1))
  df$number_of_neighbours <- lengths(adj)
  df
}


#' Estimate Facet Size from Neighbour Centre Spacing
#'
#' Estimates local facet size from the three-dimensional centre-to-centre
#' spacing between neighbouring facets. CV3D currently assumes all coordinate
#' columns are expressed in micrometres (µm).
#'
#' For each facet, the mean Euclidean centre-to-centre distance to all available
#' neighbours is used as the initial facet-size estimate. This is a geometric
#' proxy for local facet diameter rather than a direct measurement of the facet
#' boundary. A spatially smoothed
#' estimate is additionally calculated by giving equal weight to the focal
#' facet's estimate and the mean size estimate of its neighbouring facets.
#'
#' @param df A data frame or tibble containing facet identifiers,
#'   three-dimensional facet-centre coordinates, and a column listing the
#'   neighbouring facet identifiers.
#' @param id_col Character. Name of the column containing unique facet
#'   identifiers. Default: `"ID"`.
#' @param x_col Character. Name of the column containing x coordinates.
#'   Default: `"x"`.
#' @param y_col Character. Name of the column containing y coordinates.
#'   Default: `"y"`.
#' @param z_col Character. Name of the column containing z coordinates.
#'   Default: `"z"`.
#' @param neighbours_col Character. Name of the column containing neighbouring
#'   facet identifiers. Default: `"neighbours"`.
#' @param sep Character. Separator used between facet identifiers in
#'   `neighbours_col`. Default: `";"`.
#' @param keep_distances Logical. If `TRUE`, return both the facet-size summary
#'   and the individual focal-to-neighbour distances. If `FALSE`, return only
#'   the summary. Default: `FALSE`.
#'
#' @return If `keep_distances = FALSE`, a data frame containing `ID`,
#'   `facet_size`, the mean Euclidean centre spacing to neighbouring facets in
#'   micrometres; `n_neighbours_used`, the number of valid neighbour distances
#'   used; and `facet_size_smoothed`, a locally smoothed estimate giving equal
#'   weight to the focal estimate and the mean estimate of its neighbours. If
#'   `keep_distances = TRUE`, a list with `summary` and `distances`.
#'
#' @examples
#' data(cv3d_example_facets)
#' neighbours <- find_neighbours(cv3d_example_facets)
#' facet_sizes <- calculate_facet_size(neighbours)
#' head(facet_sizes)
#'
#' @export
calculate_facet_size <- function(df,
                                 id_col = "ID",
                                 x_col = "x",
                                 y_col = "y",
                                 z_col = "z",
                                 neighbours_col = "neighbours",
                                 sep = ";",
                                 keep_distances = FALSE) {
  required <- c(id_col, x_col, y_col, z_col, neighbours_col)
  if (!all(required %in% names(df))) stop("Missing required columns.", call. = FALSE)
  if (anyDuplicated(df[[id_col]]) > 0) stop("IDs must be unique.", call. = FALSE)

  ids <- as.character(df[[id_col]])
  coords <- as.matrix(df[, c(x_col, y_col, z_col), drop = FALSE])
  storage.mode(coords) <- "double"

  parse_neighbours <- function(value) {
    if (length(value) == 0 || is.na(value) || !nzchar(trimws(as.character(value)))) return(character(0))
    out <- trimws(strsplit(as.character(value), split = sep, fixed = TRUE)[[1]])
    out[nzchar(out)]
  }
  neighbours <- lapply(df[[neighbours_col]], parse_neighbours)

  distance_rows <- vector("list", nrow(df))
  facet_size <- rep(NA_real_, nrow(df))
  n_used <- integer(nrow(df))

  for (i in seq_len(nrow(df))) {
    nb <- neighbours[[i]]
    idx <- match(nb, ids)
    valid <- !is.na(idx)
    idx <- idx[valid]
    nb_valid <- nb[valid]

    if (length(idx) == 0) {
      distance_rows[[i]] <- data.frame(
        ID = character(), neighbour_ID = character(),
        x = numeric(), y = numeric(), z = numeric(),
        x_nbr = numeric(), y_nbr = numeric(), z_nbr = numeric(),
        dist = numeric(), stringsAsFactors = FALSE
      )
      next
    }

    diffs <- sweep(coords[idx, , drop = FALSE], 2, coords[i, ], "-")
    distances <- sqrt(rowSums(diffs^2))
    finite <- is.finite(distances)
    facet_size[i] <- if (any(finite)) mean(distances[finite]) else NA_real_
    n_used[i] <- sum(finite)

    distance_rows[[i]] <- data.frame(
      ID = rep(ids[i], length(idx)),
      neighbour_ID = nb_valid,
      x = rep(coords[i, 1], length(idx)),
      y = rep(coords[i, 2], length(idx)),
      z = rep(coords[i, 3], length(idx)),
      x_nbr = coords[idx, 1],
      y_nbr = coords[idx, 2],
      z_nbr = coords[idx, 3],
      dist = distances,
      stringsAsFactors = FALSE
    )
  }

  neighbour_mean <- rep(NA_real_, nrow(df))
  for (i in seq_len(nrow(df))) {
    idx <- match(neighbours[[i]], ids)
    vals <- facet_size[idx[!is.na(idx)]]
    vals <- vals[is.finite(vals)]
    if (length(vals) > 0) neighbour_mean[i] <- mean(vals)
  }

  out <- data.frame(
    ID = ids,
    facet_size = facet_size,
    n_neighbours_used = n_used,
    facet_size_smoothed = (facet_size + neighbour_mean) / 2,
    stringsAsFactors = FALSE
  )

  if (keep_distances) {
    dists <- if (length(distance_rows) > 0) do.call(rbind, distance_rows) else data.frame()
    return(list(summary = out, distances = dists))
  }
  out
}


#' Estimate Facet Surface Normals
#'
#' Estimates an outward-pointing surface normal for each compound-eye facet
#' from the three-dimensional positions of the focal facet and its neighbouring
#' facets.
#'
#' For each focal facet, local triangles are formed from the focal facet centre
#' and pairs of neighbouring facet centres. Unit normals are calculated for
#' these triangles and oriented outward using the vector from the centroid of
#' the eye point cloud towards the focal facet as a reference. The local triangle
#' normals are averaged to obtain the focal estimate. The focal normal is then
#' averaged with the normals of valid neighbouring facets to reduce local
#' directional noise, and the resulting vector is renormalised to unit length.
#'
#' Facets with fewer than two valid neighbours cannot define a local surface
#' normal and are removed. Because their removal can reduce the neighbour count
#' of other facets, this filtering is repeated until all retained facets have
#' at least two retained neighbours.
#'
#' @param df A data frame or tibble containing one row per facet. It must
#'   contain facet identifiers (`ID`), facet-centre coordinates (`x`, `y`,
#'   `z`), a `neighbours` column containing neighbouring facet identifiers, and
#'   `number_of_neighbours`.
#' @param cores Integer. Number of processor cores used for parallel normal
#'   calculation. Default: `1`.
#' @param plot_file Character or `NULL`. If supplied, write a PDF containing
#'   histograms of the x, y, and z components of the calculated normals.
#'   Default: `NULL`.
#' @param plot_results Logical. If `TRUE`, plot histograms of the normal-vector
#'   components to the active graphics device. Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress and timing information.
#'   Default: `FALSE`.
#'
#' @return A tibble with one row per retained facet and the columns `ID`,
#'   `norm.x`, `norm.y`, and `norm.z`, describing the estimated outward facet
#'   normal.
#'
#' @examples
#' data(cv3d_example_facets)
#' neighbours <- find_neighbours(cv3d_example_facets)
#' facet_normals <- get_facet_normals(neighbours, cores = 1, verbose = FALSE)
#' head(facet_normals)
#'
#' @export
get_facet_normals <- function(df,
                              cores = 1,
                              plot_file = NULL,
                              plot_results = FALSE,
                              verbose = FALSE) {
  required <- c("ID", "x", "y", "z", "neighbours", "number_of_neighbours")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(cores) || length(cores) != 1 || !is.finite(cores) || cores < 1) {
    stop("cores must be a positive integer.", call. = FALSE)
  }
  cores <- as.integer(cores)

  parse_neighbours <- function(value) {
    if (is.na(value) || !nzchar(trimws(as.character(value)))) return(character(0))
    out <- trimws(strsplit(as.character(value), split = ";", fixed = TRUE)[[1]])
    out[nzchar(out)]
  }

  work <- as.data.frame(df, stringsAsFactors = FALSE)
  work$ID <- as.character(work$ID)

  repeat {
    keep_ids <- work$ID
    parsed <- lapply(work$neighbours, parse_neighbours)
    parsed <- lapply(parsed, function(v) v[v %in% keep_ids])
    work$neighbours <- vapply(parsed, paste, collapse = "; ", FUN.VALUE = character(1))
    work$number_of_neighbours <- lengths(parsed)
    keep <- work$number_of_neighbours >= 2
    if (all(keep)) break
    work <- work[keep, , drop = FALSE]
    if (nrow(work) == 0) break
  }

  if (nrow(work) == 0) {
    return(tibble::tibble(ID = character(), norm.x = numeric(), norm.y = numeric(), norm.z = numeric()))
  }

  eye_mean <- colMeans(as.matrix(work[, c("x", "y", "z"), drop = FALSE]))
  ids <- work$ID
  coords <- as.matrix(work[, c("x", "y", "z"), drop = FALSE])
  storage.mode(coords) <- "double"

  calculate_one <- function(i) {
    focal <- coords[i, ]
    centre_vector <- focal - eye_mean
    centre_length <- sqrt(sum(centre_vector^2))
    if (!is.finite(centre_length) || centre_length == 0) {
      return(c(norm.x = NA_real_, norm.y = NA_real_, norm.z = NA_real_))
    }
    centre_vector <- centre_vector / centre_length

    neighbour_ids <- parse_neighbours(work$neighbours[i])
    neighbour_idx <- match(neighbour_ids, ids)
    neighbour_idx <- neighbour_idx[!is.na(neighbour_idx)]
    if (length(neighbour_idx) < 2) {
      return(c(norm.x = NA_real_, norm.y = NA_real_, norm.z = NA_real_))
    }

    # Construct each unordered neighbour pair exactly once.  The shortest
    # neighbour-neighbour pairs approximate the polygon edges around the focal
    # facet; equal distances must remain distinct because regular lattices
    # legitimately contain several different pairs with the same length.
    pair_matrix <- utils::combn(neighbour_idx, 2)
    pairs <- data.frame(n1 = pair_matrix[1, ], n2 = pair_matrix[2, ])
    pair_dist <- sqrt(rowSums((coords[pairs$n2, , drop = FALSE] -
                               coords[pairs$n1, , drop = FALSE])^2))
    valid_pairs <- is.finite(pair_dist) & pair_dist > 0
    pairs <- pairs[valid_pairs, , drop = FALSE]
    pair_dist <- pair_dist[valid_pairs]
    if (length(pair_dist) == 0) {
      return(c(norm.x = NA_real_, norm.y = NA_real_, norm.z = NA_real_))
    }

    order_idx <- order(pair_dist)
    pairs <- pairs[order_idx, , drop = FALSE]
    pairs <- pairs[seq_len(min(nrow(pairs), length(neighbour_idx))), , drop = FALSE]

    local_normals <- matrix(NA_real_, nrow = nrow(pairs), ncol = 3)
    for (j in seq_len(nrow(pairs))) {
      curr_normal <- calculate_normal(
        A = focal,
        B = coords[pairs$n1[j], ],
        C = coords[pairs$n2[j], ],
        normalize = TRUE
      )
      curr_length <- sqrt(sum(curr_normal^2))
      if (!all(is.finite(curr_normal)) || !is.finite(curr_length) || curr_length == 0) next
      if (sum(curr_normal * centre_vector) < 0) curr_normal <- -curr_normal
      local_normals[j, ] <- curr_normal
    }

    # Average all valid local triangle normals and renormalise the focal
    # estimate. This gives each valid focal estimate equal weight during the
    # subsequent neighbour smoothing, irrespective of local triangle spread.
    valid_local <- stats::complete.cases(local_normals) & rowSums(local_normals^2) > 0
    if (!any(valid_local)) {
      return(c(norm.x = NA_real_, norm.y = NA_real_, norm.z = NA_real_))
    }
    focal_normal <- colMeans(local_normals[valid_local, , drop = FALSE])
    focal_length <- sqrt(sum(focal_normal^2))
    if (!all(is.finite(focal_normal)) || !is.finite(focal_length) || focal_length == 0) {
      return(c(norm.x = NA_real_, norm.y = NA_real_, norm.z = NA_real_))
    }
    focal_normal / focal_length
  }

  if (verbose) message("Calculating ", nrow(work), " facet normals.")

  if (cores > 1 && nrow(work) > 1) {
    cores <- min(cores, nrow(work))
    doParallel::registerDoParallel(cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    raw_list <- foreach::foreach(i = seq_len(nrow(work)), .packages = "CV3D") %dopar% {
      calculate_one(i)
    }
  } else {
    raw_list <- lapply(seq_len(nrow(work)), calculate_one)
  }

  raw_normals <- do.call(rbind, raw_list)
  colnames(raw_normals) <- c("norm.x", "norm.y", "norm.z")

  smoothed <- matrix(NA_real_, nrow = nrow(work), ncol = 3)
  for (i in seq_len(nrow(work))) {
    neighbour_idx <- match(parse_neighbours(work$neighbours[i]), ids)
    neighbour_idx <- neighbour_idx[!is.na(neighbour_idx)]
    values <- raw_normals[c(i, neighbour_idx), , drop = FALSE]
    mean_normal <- colMeans(values, na.rm = TRUE)
    normal_length <- sqrt(sum(mean_normal^2))
    if (all(is.finite(mean_normal)) && is.finite(normal_length) && normal_length > 0) {
      smoothed[i, ] <- mean_normal / normal_length
    }
  }

  result <- tibble::tibble(
    ID = ids,
    norm.x = smoothed[, 1],
    norm.y = smoothed[, 2],
    norm.z = smoothed[, 3]
  )

  plot_histograms <- function() {
    graphics::par(mfrow = c(3, 1))
    on.exit(graphics::par(mfrow = c(1, 1)), add = TRUE)
    for (component in c("norm.x", "norm.y", "norm.z")) {
      values <- result[[component]]
      values <- values[is.finite(values)]
      if (length(values) > 0) {
        graphics::hist(values, breaks = 15, main = component, xlab = component)
      }
    }
  }

  if (plot_results) plot_histograms()
  if (!is.null(plot_file)) {
    grDevices::pdf(plot_file, onefile = TRUE, paper = "a4")
    plot_histograms()
    grDevices::dev.off()
  }

  if (verbose) message("Facet normals calculated.")
  result
}


#' Calculate Local Optical Properties of Compound-Eye Facets
#'
#' Calculates local optical-property estimates from facet sizes, facet surface
#' normals, and neighbourhood relationships.
#'
#' For each facet, the inter-facet angle is the mean angular separation between
#' its surface normal and the surface normals of its valid immediate neighbours.
#' The sampling lattice is then used for Snyder's lattice-dependent sampling
#' frequency and eye parameter. For a square lattice,
#' \eqn{v_s = 1 / (2 \Delta\phi)} and \eqn{p = D\Delta\phi}. For a hexagonal
#' lattice, \eqn{v_s = 1 / (\sqrt{3}\Delta\phi)} and
#' \eqn{p = D\sqrt{3}\Delta\phi/2}. Here, \eqn{D} is the smoothed facet-size
#' estimate in micrometres and \eqn{\Delta\phi} is in radians; the resulting
#' eye parameter is therefore expressed in µm·rad.
#'
#' `acuity_cpd` is calculated separately using the interommatidial-angle
#' relationship used by Feller et al. (2020), \eqn{1 / (2\Delta\phi)}, with
#' \eqn{\Delta\phi} expressed in degrees. Here it is a local anatomical
#' acuity estimate for each facet based on that facet's mean immediate-neighbour
#' inter-facet angle; it is independent of the square/hexagonal lattice option.
#'
#' @param df A data frame or tibble containing one row per facet. It must
#'   contain facet identifiers (`ID`), neighbouring facet identifiers
#'   (`neighbours`), facet surface normals (`norm.x`, `norm.y`, `norm.z`), and
#'   smoothed facet-size estimates (`facet_size_smoothed`) in micrometres.
#' @param lattice Character. Local sampling lattice used for Snyder's sampling
#'   frequency and eye parameter. Either `"hexagonal"` (default) or `"square"`.
#' @param cores Integer. Number of processor cores used for parallel
#'   calculation. Default: `1`.
#' @param plot_results Logical. If `TRUE`, plot diagnostic histograms of the
#'   calculated optical properties. Default: `FALSE`.
#' @param plot_file Character or `NULL`. If supplied, write diagnostic
#'   histograms to a PDF file. Default: `NULL`.
#' @param verbose Logical. If `TRUE`, print progress and timing information.
#'   Default: `FALSE`.
#'
#' @return A tibble with one row per facet and the columns `ID`,
#'   `interfacet_angle_deg`, the mean immediate-neighbour inter-facet angle in
#'   degrees; `interfacet_angle_rad`, the same angle in radians;
#'   `sampling_lattice`, the selected lattice; `eye_parameter`, Snyder's local
#'   lattice-corrected eye parameter in µm·rad; `sampling_frequency_rad`,
#'   Snyder's lattice-corrected sampling frequency in cycles per radian; and
#'   `acuity_cpd`, the local anatomical acuity estimate in cycles per degree
#'   calculated from interommatidial angle using the Feller et al. relationship.
#'
#' @references
#' Snyder, A. W. (1977). Acuity of compound eyes: Physical limitations and
#' design. Journal of Comparative Physiology A, 116, 161-182.
#' \doi{10.1007/BF00605401}
#'
#' Snyder, A. W., Stavenga, D. G., & Laughlin, S. B. (1977). Spatial
#' information capacity of compound eyes. Journal of Comparative Physiology A,
#' 116, 183-207. \doi{10.1007/BF00605402}
#'
#' Feller, K. D., Sharkey, C. R., McDuffee-Altekruse, A., et al. (2020).
#' Surf and turf vision: Patterns and predictors of visual acuity in compound
#' eye evolution. Arthropod Structure & Development, 60, 101002.
#' \doi{10.1016/j.asd.2020.101002}
#'
#' @examples
#' data(cv3d_example_facets)
#' facets <- find_neighbours(cv3d_example_facets)
#' sizes <- calculate_facet_size(facets)
#' facets <- merge(facets, sizes, by = "ID", all.x = TRUE, sort = FALSE)
#' normals <- get_facet_normals(facets, cores = 1, verbose = FALSE)
#' facets <- merge(facets, normals, by = "ID", all.x = TRUE, sort = FALSE)
#' optics <- get_optic_properties(
#'   facets,
#'   lattice = "hexagonal",
#'   cores = 1,
#'   verbose = FALSE
#' )
#' head(optics)
#'
#' @export
get_optic_properties <- function(df,
                                 lattice = c("hexagonal", "square"),
                                 cores = 1,
                                 plot_results = FALSE,
                                 plot_file = NULL,
                                 verbose = FALSE) {
  required_cols <- c("ID", "neighbours", "norm.x", "norm.y", "norm.z", "facet_size_smoothed")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  lattice <- match.arg(lattice)
  if (!is.numeric(cores) || length(cores) != 1 || !is.finite(cores) || cores < 1) {
    stop("cores must be a positive integer.", call. = FALSE)
  }
  cores <- as.integer(cores)

  ids <- as.character(df$ID)
  parse_neighbours <- function(value) {
    if (is.na(value) || !nzchar(trimws(as.character(value)))) return(character(0))
    out <- trimws(strsplit(as.character(value), split = ";", fixed = TRUE)[[1]])
    out[nzchar(out)]
  }

  calculate_one <- function(i) {
    focal_normal <- as.numeric(df[i, c("norm.x", "norm.y", "norm.z")])
    focal_valid <- all(is.finite(focal_normal)) && sqrt(sum(focal_normal^2)) > 0
    angles <- numeric(0)

    if (focal_valid) {
      neighbour_ids <- parse_neighbours(df$neighbours[i])
      neighbour_idx <- match(neighbour_ids, ids)
      neighbour_idx <- neighbour_idx[!is.na(neighbour_idx)]
      for (idx in neighbour_idx) {
        neighbour_normal <- as.numeric(df[idx, c("norm.x", "norm.y", "norm.z")])
        if (!all(is.finite(neighbour_normal)) || sqrt(sum(neighbour_normal^2)) == 0) next
        angle <- vector_angle(focal_normal, neighbour_normal, unit = "radians")
        if (is.finite(angle)) angles <- c(angles, angle)
      }
    }

    if (length(angles) > 0) mean(angles) else NA_real_
  }

  if (verbose) message("Calculating inter-facet angles for ", nrow(df), " facets.")

  if (cores > 1 && nrow(df) > 1) {
    cores <- min(cores, nrow(df))
    doParallel::registerDoParallel(cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)
    i <- NULL  # foreach iteration variable; explicit binding for R CMD check
    angle_list <- foreach::foreach(i = seq_len(nrow(df)), .packages = "CV3D") %dopar% {
      calculate_one(i)
    }
  } else {
    angle_list <- lapply(seq_len(nrow(df)), calculate_one)
  }

  interfacet_angle_rad <- as.numeric(unlist(angle_list, use.names = FALSE))
  interfacet_angle_deg <- interfacet_angle_rad * 180 / pi
  facet_size_smoothed <- as.numeric(df$facet_size_smoothed)

  if (lattice == "hexagonal") {
    eye_parameter <- facet_size_smoothed * interfacet_angle_rad * (sqrt(3) / 2)
    sampling_frequency_rad <- 1 / (sqrt(3) * interfacet_angle_rad)
  } else {
    eye_parameter <- facet_size_smoothed * interfacet_angle_rad
    sampling_frequency_rad <- 1 / (2 * interfacet_angle_rad)
  }

  # Local anatomical acuity estimate using the interommatidial-angle
  # relationship used by Feller et al. (2020).
  acuity_cpd <- 1 / (2 * interfacet_angle_deg)

  invalid_angle <- !is.finite(interfacet_angle_rad) | interfacet_angle_rad <= 0
  invalid_size <- !is.finite(facet_size_smoothed) | facet_size_smoothed <= 0

  eye_parameter[invalid_angle | invalid_size] <- NA_real_
  sampling_frequency_rad[invalid_angle] <- NA_real_
  acuity_cpd[invalid_angle] <- NA_real_

  result <- tibble::tibble(
    ID = ids,
    interfacet_angle_deg = interfacet_angle_deg,
    interfacet_angle_rad = interfacet_angle_rad,
    sampling_lattice = rep(lattice, length(ids)),
    eye_parameter = eye_parameter,
    sampling_frequency_rad = sampling_frequency_rad,
    acuity_cpd = acuity_cpd
  )

  plot_hist <- function(x, main, xlab) {
    values <- x[is.finite(x)]
    if (length(values) > 0) graphics::hist(values, breaks = 15, main = main, xlab = xlab)
  }
  plot_all <- function() {
    graphics::par(mfrow = c(5, 1))
    on.exit(graphics::par(mfrow = c(1, 1)), add = TRUE)
    plot_hist(facet_size_smoothed, "Facet diameter", "Facet diameter")
    plot_hist(result$interfacet_angle_deg, "Inter-facet angles", "Inter-facet angle (degrees)")
    plot_hist(result$eye_parameter, "Eye parameter", paste0("Eye parameter (", lattice, " lattice)"))
    plot_hist(result$sampling_frequency_rad, "Sampling frequency", "Cycles per radian")
    plot_hist(result$acuity_cpd, "Local anatomical acuity", "Cycles per degree")
  }

  if (plot_results) plot_all()
  if (!is.null(plot_file)) {
    grDevices::pdf(plot_file, onefile = TRUE, paper = "a4")
    plot_all()
    grDevices::dev.off()
  }

  if (verbose) message("Optical properties calculated.")
  result
}


#' Convert Cartesian Coordinates to View Angles
#'
#' Converts three-dimensional Cartesian coordinates to elevation and azimuth
#' relative to a specified centre.
#'
#' Coordinates are first translated so that `center` becomes the origin.
#' Elevation is measured from the x-y plane, with positive values towards the
#' positive z-axis. By default, azimuth follows the CV3D anatomical convention
#' in which the anterior viewing direction (negative y-axis after landmark
#' alignment) is 0 degrees, the positive x-axis is +90 degrees, and the
#' negative x-axis is -90 degrees.
#'
#' @param x Numeric vector of x coordinates.
#' @param y Numeric vector of y coordinates.
#' @param z Numeric vector of z coordinates.
#' @param center Numeric vector of length three specifying the x, y, and z
#'   coordinates of the angular reference centre. Default: `c(0, 0, 0)`.
#' @param front_zero Logical. If `TRUE`, rotate azimuth by 90 degrees so that
#'   the negative y-axis, corresponding to the anterior direction in the CV3D
#'   global coordinate system, has azimuth 0 degrees. If `FALSE`, use the
#'   conventional Cartesian azimuth measured from the positive x-axis.
#'   Default: `TRUE`.
#'
#' @return A tibble with the columns `elevation` and `azimuth`, both in degrees.
#'   Azimuth is wrapped to the interval from -180 degrees inclusive to +180
#'   degrees exclusive. Coordinates located exactly at `center` return `NA`
#'   because their angular direction is undefined.
#'
#' @export
#' @examples
#' directions <- cartesian_to_view_angles(
#'   x = c(0, 1, 0, -1),
#'   y = c(-1, 0, 1, 0),
#'   z = c(0, 0, 0, 0)
#' )
#' directions
cartesian_to_view_angles <- function(x,
                                     y,
                                     z,
                                     center = c(0, 0, 0),
                                     front_zero = TRUE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  z <- as.numeric(z)
  center <- as.numeric(center)

  if (length(x) != length(y) || length(x) != length(z)) {
    stop("x, y, and z must have equal lengths.", call. = FALSE)
  }
  if (length(center) != 3 || any(!is.finite(center))) {
    stop("center must contain three finite numeric values.", call. = FALSE)
  }
  if (length(front_zero) != 1 || is.na(front_zero) || !is.logical(front_zero)) {
    stop("front_zero must be TRUE or FALSE.", call. = FALSE)
  }

  dx <- x - center[1]
  dy <- y - center[2]
  dz <- z - center[3]
  radius <- sqrt(dx^2 + dy^2 + dz^2)

  elevation <- rep(NA_real_, length(radius))
  azimuth <- rep(NA_real_, length(radius))
  valid <- is.finite(dx) & is.finite(dy) & is.finite(dz) & is.finite(radius) & radius > 0

  elevation[valid] <- asin(pmax(-1, pmin(1, dz[valid] / radius[valid]))) * 180 / pi
  azimuth[valid] <- atan2(dy[valid], dx[valid]) * 180 / pi
  if (front_zero) azimuth[valid] <- azimuth[valid] + 90
  azimuth[valid] <- ((azimuth[valid] + 180) %% 360) - 180

  tibble::tibble(elevation = elevation, azimuth = azimuth)
}


#' Calculate the Intersection of a 3D Ray and Sphere
#'
#' Calculates the nearest forward intersection between a ray in
#' three-dimensional space and a sphere. The ray originates at `point` and
#' extends in `direction`.
#'
#' In CV3D, this function is used to project facet surface normals from facet
#' centres onto a corneal-projection sphere.
#'
#' @param point Numeric vector of length three containing the x, y, and z
#'   coordinates of the ray origin.
#' @param direction Numeric vector of length three giving the x, y, and z
#'   components of the ray direction. The vector does not need to have unit
#'   length.
#' @param sphere_center Numeric vector of length three containing the x, y, and
#'   z coordinates of the sphere centre.
#' @param sphere_radius Numeric. Radius of the sphere. Must be greater than zero
#'   and use the same spatial units as `point` and `sphere_center`.
#'
#' @return A numeric vector of length three containing the x, y, and z
#'   coordinates of the nearest sphere intersection in the forward ray
#'   direction. If no forward intersection exists, all three values are
#'   `NA_real_`.
#'
#' @export
#' @examples
#' ray_sphere_intersection(
#'   point = c(0, 0, 0),
#'   direction = c(1, 0, 0),
#'   sphere_center = c(0, 0, 0),
#'   sphere_radius = 10
#' )
ray_sphere_intersection <- function(point,
                                    direction,
                                    sphere_center,
                                    sphere_radius) {
  point <- as.numeric(point)
  direction <- as.numeric(direction)
  sphere_center <- as.numeric(sphere_center)
  sphere_radius <- as.numeric(sphere_radius)

  if (length(point) != 3 || any(!is.finite(point))) {
    stop("point must contain three finite numeric values.", call. = FALSE)
  }
  if (length(direction) != 3 || any(!is.finite(direction)) || sqrt(sum(direction^2)) == 0) {
    stop("direction must be a non-zero vector containing three finite numeric values.", call. = FALSE)
  }
  if (length(sphere_center) != 3 || any(!is.finite(sphere_center))) {
    stop("sphere_center must contain three finite numeric values.", call. = FALSE)
  }
  if (length(sphere_radius) != 1 || !is.finite(sphere_radius) || sphere_radius <= 0) {
    stop("sphere_radius must be a single positive finite number.", call. = FALSE)
  }

  offset <- point - sphere_center
  A <- sum(direction^2)
  B <- 2 * sum(offset * direction)
  C <- sum(offset^2) - sphere_radius^2
  discriminant <- B^2 - 4 * A * C

  tolerance <- 100 * .Machine$double.eps * max(1, abs(B^2), abs(4 * A * C))
  if (discriminant < -tolerance) {
    return(c(x = NA_real_, y = NA_real_, z = NA_real_))
  }
  if (discriminant < 0) discriminant <- 0

  root <- sqrt(discriminant)
  solutions <- c((-B - root) / (2 * A), (-B + root) / (2 * A))
  forward <- solutions[is.finite(solutions) & solutions > 0]
  if (length(forward) == 0) {
    return(c(x = NA_real_, y = NA_real_, z = NA_real_))
  }

  t <- min(forward)
  result <- point + t * direction
  stats::setNames(as.numeric(result), c("x", "y", "z"))
}
