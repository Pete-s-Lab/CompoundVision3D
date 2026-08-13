#' Detect Compound-Eye Edge Facets from Local Angular Coverage
#'
#' Detects facets at the boundary of a sampled compound-eye surface from the
#' largest empty angular sector around each facet in a local tangent plane.
#' Coordinates are assumed to be in micrometres (µm), consistent with CV3D.
#'
#' @param df Data frame or tibble containing facet IDs and x/y/z coordinates.
#' @param gap_threshold_deg Numeric. A facet is classified as an edge facet when
#'   its largest local angular gap is greater than this threshold, in degrees.
#'   Default: `90`.
#' @param x,y,z,id Character. Column names for coordinates and facet IDs.
#' @param k Positive integer. Number of nearest facet centres used to assess
#'   angular coverage. Default: `6`.
#' @param tangent_k Positive integer. Number of nearby facets used to estimate
#'   the local tangent plane. Default: `12`.
#' @param knn_search Positive integer. Number of nearest neighbours retained for
#'   local calculations. Default: `20`.
#'
#' @details
#' This function performs edge classification only; it does not assign or
#' remove facet-neighbour links. For each focal facet, CV3D first estimates a
#' local tangent plane from the focal facet and its nearby facet centres. The
#' `k` closest candidate centres are projected into that plane and converted to
#' polar angles. After sorting those angles around 360 degrees, the largest
#' circular gap - including the wrap-around gap - is stored as
#' `edge_angular_gap_deg`.
#'
#' A well-surrounded interior facet normally has candidate directions covering
#' the full circle, whereas a boundary facet has a conspicuous unsampled
#' sector. A facet is therefore classified as an edge when its largest gap is
#' strictly greater than `gap_threshold_deg`. The threshold is deliberately an
#' explicit QC parameter rather than an automatically fitted biological
#' boundary. The CV3D UI can compare several thresholds before the final
#' neighbour graph is generated.
#'
#' The tangent plane is obtained by singular-value decomposition of the local
#' centred coordinate cloud. Consequently, edge detection is based on local 3D
#' geometry rather than on any particular global eye orientation. Translation
#' and rigid rotation of the eye do not change the angular-gap criterion.
#' Coordinates are assumed to be in micrometres (µm), as elsewhere in CV3D,
#' although the reported angular gap itself is dimensionless and expressed in
#' degrees.
#'
#' @return A data frame with one row per input facet and columns `ID`,
#'   `edge_angular_gap_deg`, and `is_edge_facet`. Input row order is preserved.
#'
#' @examples
#' data(cv3d_example_facets)
#' edge_qc <- detect_facet_edges(cv3d_example_facets, gap_threshold_deg = 90)
#' table(edge_qc$is_edge_facet)
#'
#' @export
detect_facet_edges <- function(df,
                               gap_threshold_deg = 90,
                               x = "x", y = "y", z = "z", id = "ID",
                               k = 6,
                               tangent_k = 12,
                               knn_search = 20) {
  prep <- .cv3d_prepare_edge_geometry(
    df = df, x = x, y = y, z = z, id = id,
    k = k, tangent_k = tangent_k, knn_search = knn_search
  )

  if (!is.numeric(gap_threshold_deg) || length(gap_threshold_deg) != 1L ||
      !is.finite(gap_threshold_deg) || gap_threshold_deg <= 0 || gap_threshold_deg >= 360) {
    stop("gap_threshold_deg must be a single finite value between 0 and 360 degrees.", call. = FALSE)
  }

  data.frame(
    ID = as.character(df[[id]]),
    edge_angular_gap_deg = prep$edge_angular_gap_deg,
    is_edge_facet = is.finite(prep$edge_angular_gap_deg) & prep$edge_angular_gap_deg > gap_threshold_deg,
    stringsAsFactors = FALSE
  )
}


#' Find Facet Neighbours with Edge-Aware Local Geometry
#'
#' Finds reciprocal first-ring facet neighbours while treating the outer eye
#' boundary explicitly. The method uses reciprocal six-nearest-neighbour
#' membership, a robust local core-spacing gate, edge detection from a large
#' empty angular sector in the local tangent plane, and conservative
#' angle-shadow pruning for already detected edge facets.
#'
#' @param df Data frame or tibble containing facet IDs and x/y/z coordinates.
#' @param edge_gap_threshold_deg Numeric edge-detection threshold in degrees.
#'   Default: `90`.
#' @param x,y,z,id Character. Column names for coordinates and facet IDs.
#' @param k Positive integer. Maximum local nearest-neighbour candidate count.
#'   Default: `6`.
#' @param core_k Positive integer. Number of closest facet-centre distances used
#'   for the robust local core spacing. Default: `3`.
#' @param tangent_k Positive integer. Number of nearby facets used to estimate
#'   the local tangent plane. Default: `12`.
#' @param knn_search Positive integer. Number of nearest neighbours retained for
#'   local calculations. Default: `20`.
#' @param interior_link_factor Positive numeric. Maximum link length relative to
#'   local core spacing when both endpoints are interior. Default: `1.5`.
#' @param edge_link_factor Positive numeric. Maximum link length relative to
#'   local core spacing when either endpoint is an edge facet. Default: `1.4`.
#' @param shadow_angle_fraction Positive numeric. Fraction of the median angular
#'   spacing in six-neighbour interior facets used for the edge shadow-angle
#'   cutoff. Default: `2/3`.
#' @param shadow_angle_min_deg,shadow_angle_max_deg Numeric bounds for the
#'   derived shadow-angle cutoff. Defaults: `30` and `45` degrees.
#' @param shadow_radial_ratio Numeric. The farther of two near-collinear edge
#'   candidates must be at least this many times farther away before it can be
#'   removed. Default: `1.15`.
#' @param shadow_min_remaining Non-negative integer. Never prune an edge facet
#'   below this number of retained neighbours. Default: `2`.
#' @param verbose Logical. Print a short method summary. Default: `FALSE`.
#'
#' @details
#' The method is designed for an approximately first-ring compound-eye lattice
#' in which interior facets may have up to `k` direct neighbours but facets at
#' the sampled eye boundary should not be filled artificially with second-ring
#' neighbours. It proceeds in the following stages:
#'
#' 1. Local tangent geometry and the largest angular gap are calculated as in
#'    [detect_facet_edges()]. Facets with a gap strictly greater than
#'    `edge_gap_threshold_deg` are flagged as edge facets.
#' 2. Candidate links are restricted to reciprocal `k`-nearest-neighbour
#'    membership: facet A must be among facet B's first `k` candidates and vice
#'    versa. This prevents one-sided long links.
#' 3. A robust local core spacing is calculated for every facet as the median
#'    distance to its `core_k` closest facet centres. Reciprocal links are kept
#'    only if their 3D length is no greater than `interior_link_factor` times
#'    the larger endpoint core spacing for an interior-interior link, or
#'    `edge_link_factor` times that spacing when either endpoint is an edge.
#' 4. CV3D derives an angular-shadow cutoff from the typical angular separation
#'    of retained neighbours at six-neighbour interior facets. The median
#'    reference separation is multiplied by `shadow_angle_fraction` and then
#'    constrained to `shadow_angle_min_deg`--`shadow_angle_max_deg`.
#' 5. Only at already detected edge facets, two retained neighbours that point
#'    in nearly the same local tangent-plane direction are treated as a
#'    possible first-ring/second-ring conflict. The farther candidate is pruned
#'    only when it is at least `shadow_radial_ratio` times farther away. Links
#'    are removed reciprocally, and pruning never reduces the focal facet below
#'    `shadow_min_remaining` retained neighbours.
#'
#' Edge facets are not forced to a predetermined neighbour count and there is
#' no refill stage. The resulting degree is therefore allowed to reflect the
#' available geometry naturally, with `k` acting as a maximum. The stored
#' `number_of_neighbours` column is intended both for downstream calculations
#' and for QC of the boundary solution.
#'
#' The `neighbours` column contains semicolon-separated facet IDs. The
#' `neighbour_core_spacing_um` column is in µm. Attributes on the returned data
#' frame record the derived shadow cutoff, its interior reference statistics,
#' the individual removed shadow links, and the method identifier so that the
#' exact neighbour construction can be traced downstream.
#'
#' @return The input data with additional columns `neighbours`,
#'   `number_of_neighbours`, `is_edge_facet`, `edge_angular_gap_deg`,
#'   `edge_gap_threshold_deg`, `neighbour_core_spacing_um`, and
#'   `shadow_links_removed`. The returned object has attributes
#'   `shadow_angle_threshold_deg`, `shadow_removed_links`, and
#'   `neighbour_method`.
#'
#' @examples
#' data(cv3d_example_facets)
#' nb <- find_neighbours_edge_aware(
#'   cv3d_example_facets,
#'   edge_gap_threshold_deg = 90
#' )
#' head(nb[, c("ID", "neighbours", "number_of_neighbours", "is_edge_facet")])
#'
#' @export
find_neighbours_edge_aware <- function(df,
                                       edge_gap_threshold_deg = 90,
                                       x = "x", y = "y", z = "z", id = "ID",
                                       k = 6,
                                       core_k = 3,
                                       tangent_k = 12,
                                       knn_search = 20,
                                       interior_link_factor = 1.5,
                                       edge_link_factor = 1.4,
                                       shadow_angle_fraction = 2 / 3,
                                       shadow_angle_min_deg = 30,
                                       shadow_angle_max_deg = 45,
                                       shadow_radial_ratio = 1.15,
                                       shadow_min_remaining = 2,
                                       verbose = FALSE) {
  prep <- .cv3d_prepare_edge_geometry(
    df = df, x = x, y = y, z = z, id = id,
    k = k, tangent_k = tangent_k, knn_search = knn_search
  )

  .cv3d_validate_scalar(edge_gap_threshold_deg, "edge_gap_threshold_deg", lower = 0, upper = 360, strict_lower = TRUE, strict_upper = TRUE)
  .cv3d_validate_scalar(interior_link_factor, "interior_link_factor", lower = 0, strict_lower = TRUE)
  .cv3d_validate_scalar(edge_link_factor, "edge_link_factor", lower = 0, strict_lower = TRUE)
  .cv3d_validate_scalar(shadow_angle_fraction, "shadow_angle_fraction", lower = 0, strict_lower = TRUE)
  .cv3d_validate_scalar(shadow_angle_min_deg, "shadow_angle_min_deg", lower = 0, strict_lower = TRUE)
  .cv3d_validate_scalar(shadow_angle_max_deg, "shadow_angle_max_deg", lower = shadow_angle_min_deg, strict_lower = FALSE)
  .cv3d_validate_scalar(shadow_radial_ratio, "shadow_radial_ratio", lower = 1, strict_lower = FALSE)

  core_k <- .cv3d_positive_integer(core_k, "core_k")
  shadow_min_remaining <- .cv3d_nonnegative_integer(shadow_min_remaining, "shadow_min_remaining")

  n <- nrow(prep$pts)
  k_eff <- min(.cv3d_positive_integer(k, "k"), ncol(prep$knn_idx))
  core_eff <- min(core_k, k_eff)
  if (core_eff < 1L) stop("At least one core-spacing neighbour is required.", call. = FALSE)

  edge_flags <- is.finite(prep$edge_angular_gap_deg) &
    prep$edge_angular_gap_deg > edge_gap_threshold_deg

  nn <- prep$knn_idx[, seq_len(k_eff), drop = FALSE]
  d <- prep$knn_dst[, seq_len(k_eff), drop = FALSE]
  core_scale <- apply(d[, seq_len(core_eff), drop = FALSE], 1, stats::median, na.rm = TRUE)

  adj <- replicate(n, integer(0), simplify = FALSE)
  for (i in seq_len(n)) {
    for (j in nn[i, ]) {
      if (!(i %in% nn[j, ])) next
      link_factor <- if (isTRUE(edge_flags[i]) || isTRUE(edge_flags[j])) edge_link_factor else interior_link_factor
      dij <- .cv3d_edge_distance(prep$pts, i, j)
      if (is.finite(dij) && dij <= link_factor * max(core_scale[i], core_scale[j])) {
        adj[[i]] <- unique(c(adj[[i]], j))
        adj[[j]] <- unique(c(adj[[j]], i))
      }
    }
  }
  adj <- lapply(adj, function(v) sort(unique(v)))

  shadow_cal <- .cv3d_derive_shadow_angle_threshold(
    adj = adj,
    edge_flags = edge_flags,
    prep = prep,
    fraction = shadow_angle_fraction,
    min_deg = shadow_angle_min_deg,
    max_deg = shadow_angle_max_deg
  )

  pruned <- .cv3d_prune_edge_angular_shadows(
    adj = adj,
    edge_flags = edge_flags,
    edge_gap_deg = prep$edge_angular_gap_deg,
    prep = prep,
    angle_threshold_deg = shadow_cal$threshold,
    radial_ratio = shadow_radial_ratio,
    min_remaining = shadow_min_remaining,
    ids = as.character(df[[id]])
  )
  adj <- pruned$adj

  ids <- as.character(df[[id]])
  out <- df
  out$neighbours <- vapply(adj, function(v) {
    if (length(v) == 0L) "" else paste(ids[v], collapse = "; ")
  }, character(1))
  out$number_of_neighbours <- lengths(adj)
  out$is_edge_facet <- edge_flags
  out$edge_angular_gap_deg <- prep$edge_angular_gap_deg
  out$edge_gap_threshold_deg <- rep(as.numeric(edge_gap_threshold_deg), n)
  out$neighbour_core_spacing_um <- core_scale

  removed_count <- integer(n)
  if (nrow(pruned$removed) > 0L) {
    focal_match <- match(pruned$removed$focal_ID, ids)
    focal_match <- focal_match[!is.na(focal_match)]
    if (length(focal_match) > 0L) {
      tab <- table(focal_match)
      removed_count[as.integer(names(tab))] <- as.integer(tab)
    }
  }
  out$shadow_links_removed <- removed_count

  attr(out, "shadow_angle_threshold_deg") <- shadow_cal$threshold
  attr(out, "shadow_reference_interior_median_gap_deg") <- shadow_cal$interior_median_gap
  attr(out, "shadow_reference_gap_count") <- shadow_cal$n_reference_gaps
  attr(out, "shadow_removed_links") <- pruned$removed
  attr(out, "neighbour_method") <- "edge_aware_mutual6_coregate_angle_shadow"

  if (isTRUE(verbose)) {
    message(
      "Edge-aware neighbours: ", sum(edge_flags), "/", n,
      " edge facets at gap > ", edge_gap_threshold_deg,
      " deg; shadow cutoff ", round(shadow_cal$threshold, 2),
      " deg; removed ", nrow(pruned$removed), " shadowed links."
    )
  }

  out
}


.cv3d_validate_scalar <- function(value, name, lower = -Inf, upper = Inf,
                                  strict_lower = FALSE, strict_upper = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
    stop(name, " must be a single finite numeric value.", call. = FALSE)
  }
  lower_bad <- if (strict_lower) value <= lower else value < lower
  upper_bad <- if (strict_upper) value >= upper else value > upper
  if (lower_bad || upper_bad) stop(name, " is outside its allowed range.", call. = FALSE)
  invisible(value)
}

.cv3d_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 1 || value != as.integer(value)) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

.cv3d_nonnegative_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value < 0 || value != as.integer(value)) {
    stop(name, " must be a non-negative integer.", call. = FALSE)
  }
  as.integer(value)
}

.cv3d_prepare_edge_geometry <- function(df, x, y, z, id, k, tangent_k, knn_search) {
  if (!is.data.frame(df)) stop("df must be a data frame or tibble.", call. = FALSE)
  if (!all(c(x, y, z, id) %in% names(df))) stop("Missing required coordinate or ID columns.", call. = FALSE)
  if (anyDuplicated(df[[id]]) > 0L) stop("IDs must be unique.", call. = FALSE)

  k <- .cv3d_positive_integer(k, "k")
  tangent_k <- .cv3d_positive_integer(tangent_k, "tangent_k")
  knn_search <- .cv3d_positive_integer(knn_search, "knn_search")
  if (knn_search < max(k, tangent_k)) stop("knn_search must be >= both k and tangent_k.", call. = FALSE)

  pts <- as.matrix(df[, c(x, y, z), drop = FALSE])
  storage.mode(pts) <- "double"
  if (nrow(pts) < 4L) stop("At least four facet positions are required.", call. = FALSE)
  if (any(!is.finite(pts))) stop("Facet coordinates must be finite.", call. = FALSE)

  n <- nrow(pts)
  knn_eff <- min(max(knn_search, tangent_k, k), n - 1L)
  if (knn_eff < 3L) stop("At least three nearby facets are required.", call. = FALSE)

  nn_obj <- RANN::nn2(pts, pts, k = knn_eff + 1L)
  knn_idx <- nn_obj$nn.idx[, -1, drop = FALSE]
  knn_dst <- nn_obj$nn.dists[, -1, drop = FALSE]
  tangent_eff <- min(tangent_k, ncol(knn_idx))
  edge_k_eff <- min(k, ncol(knn_idx))

  basis_list <- vector("list", n)
  max_gap <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    local_idx <- c(i, knn_idx[i, seq_len(tangent_eff)])
    local_pts <- pts[local_idx, , drop = FALSE]
    local_centred <- sweep(local_pts, 2, colMeans(local_pts), "-")
    s <- tryCatch(base::svd(local_centred), error = function(e) NULL)
    if (is.null(s) || ncol(s$v) < 2L) next
    basis <- s$v[, 1:2, drop = FALSE]
    basis_list[[i]] <- basis

    cand <- knn_idx[i, seq_len(edge_k_eff)]
    a <- .cv3d_angles_to_candidates_deg(pts, i, cand, basis)
    gaps <- .cv3d_circular_gaps_deg(a)
    if (length(gaps) > 0L) max_gap[i] <- max(gaps)
  }

  list(
    pts = pts,
    knn_idx = knn_idx,
    knn_dst = knn_dst,
    basis = basis_list,
    edge_angular_gap_deg = max_gap
  )
}

.cv3d_angles_to_candidates_deg <- function(pts, i, candidate_idx, basis) {
  candidate_idx <- as.integer(candidate_idx)
  candidate_idx <- candidate_idx[
    is.finite(candidate_idx) & candidate_idx >= 1L & candidate_idx <= nrow(pts) & candidate_idx != i
  ]
  if (length(candidate_idx) == 0L || is.null(basis)) return(numeric(0))
  vec <- pts[candidate_idx, , drop = FALSE] -
    matrix(pts[i, ], nrow = length(candidate_idx), ncol = 3L, byrow = TRUE)
  uv <- vec %*% basis
  (atan2(uv[, 2], uv[, 1]) * 180 / pi) %% 360
}

.cv3d_circular_gaps_deg <- function(angles_deg) {
  a <- sort(angles_deg[is.finite(angles_deg)])
  if (length(a) < 2L) return(numeric(0))
  diff(c(a, a[1] + 360))
}

.cv3d_edge_distance <- function(pts, i, j) {
  sqrt(sum((pts[i, ] - pts[j, ])^2))
}

.cv3d_derive_shadow_angle_threshold <- function(adj, edge_flags, prep,
                                                fraction, min_deg, max_deg) {
  interior_gaps <- numeric(0)
  for (i in seq_along(adj)) {
    if (isTRUE(edge_flags[i]) || length(adj[[i]]) != 6L) next
    a <- .cv3d_angles_to_candidates_deg(prep$pts, i, adj[[i]], prep$basis[[i]])
    g <- .cv3d_circular_gaps_deg(a)
    if (length(g) == 6L && all(is.finite(g))) interior_gaps <- c(interior_gaps, g)
  }

  if (length(interior_gaps) < 12L) {
    return(list(threshold = 40, interior_median_gap = NA_real_, n_reference_gaps = length(interior_gaps)))
  }

  med <- stats::median(interior_gaps, na.rm = TRUE)
  threshold <- fraction * med
  threshold <- max(min_deg, min(max_deg, threshold))
  list(threshold = threshold, interior_median_gap = med, n_reference_gaps = length(interior_gaps))
}

.cv3d_prune_edge_angular_shadows <- function(adj, edge_flags, edge_gap_deg, prep,
                                             angle_threshold_deg, radial_ratio,
                                             min_remaining, ids) {
  adj <- lapply(adj, function(v) sort(unique(v)))
  removed <- list()
  removal_index <- 0L

  repeat {
    changed <- FALSE
    for (i in which(edge_flags)) {
      if (length(adj[[i]]) <= min_remaining) next
      cand <- adj[[i]]
      if (length(cand) < 2L) next

      a <- .cv3d_angles_to_candidates_deg(prep$pts, i, cand, prep$basis[[i]])
      if (length(a) != length(cand) || any(!is.finite(a))) next
      ord <- order(a)
      a <- a[ord]
      cand <- cand[ord]
      gaps <- diff(c(a, a[1] + 360))

      choices <- list()
      choice_index <- 0L
      for (kk in seq_along(gaps)) {
        gap <- gaps[kk]
        if (!is.finite(gap) || gap >= angle_threshold_deg) next

        j1 <- cand[kk]
        j2 <- cand[if (kk == length(cand)) 1L else kk + 1L]
        d1 <- .cv3d_edge_distance(prep$pts, i, j1)
        d2 <- .cv3d_edge_distance(prep$pts, i, j2)
        near_d <- min(d1, d2)
        far_d <- max(d1, d2)
        if (!is.finite(near_d) || near_d <= 0 || !is.finite(far_d)) next
        ratio <- far_d / near_d
        if (ratio < radial_ratio) next

        drop <- if (d1 > d2) j1 else j2
        keep <- if (drop == j1) j2 else j1
        choice_index <- choice_index + 1L
        choices[[choice_index]] <- list(
          gap = gap, ratio = ratio, drop = drop, keep = keep,
          near_distance = near_d, far_distance = far_d
        )
      }
      if (length(choices) == 0L) next

      choice_gaps <- vapply(choices, `[[`, numeric(1), "gap")
      choice_ratios <- vapply(choices, `[[`, numeric(1), "ratio")
      z <- choices[[order(choice_gaps, -choice_ratios)[1]]]
      if (length(adj[[i]]) <= min_remaining || !(z$drop %in% adj[[i]])) next

      adj[[i]] <- setdiff(adj[[i]], z$drop)
      adj[[z$drop]] <- setdiff(adj[[z$drop]], i)
      changed <- TRUE

      removal_index <- removal_index + 1L
      removed[[removal_index]] <- data.frame(
        focal_ID = ids[i],
        removed_neighbour_ID = ids[z$drop],
        retained_competing_ID = ids[z$keep],
        angular_separation_deg = z$gap,
        nearer_distance_um = z$near_distance,
        removed_distance_um = z$far_distance,
        distance_ratio = z$ratio,
        focal_edge_gap_deg = edge_gap_deg[i],
        stringsAsFactors = FALSE
      )
    }
    if (!changed) break
  }

  removed_df <- if (length(removed) > 0L) do.call(rbind, removed) else data.frame(
    focal_ID = character(0),
    removed_neighbour_ID = character(0),
    retained_competing_ID = character(0),
    angular_separation_deg = numeric(0),
    nearer_distance_um = numeric(0),
    removed_distance_um = numeric(0),
    distance_ratio = numeric(0),
    focal_edge_gap_deg = numeric(0),
    stringsAsFactors = FALSE
  )

  list(adj = lapply(adj, function(v) sort(unique(v))), removed = removed_df)
}
