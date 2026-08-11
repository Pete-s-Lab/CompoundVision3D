


#' Condense Local-Height Points into Facet Candidates
#'
#' Identifies facet candidates from a thresholded local-height point cloud
#' using iterative, height-weighted spatial condensation. At each iteration,
#' every point moves towards the height-weighted centre of neighbouring points
#' within a spherical Euclidean neighbourhood. Iteration continues until the
#' maximum displacement falls below a specified tolerance or the maximum
#' number of iterations is reached.
#'
#' After condensation, converged points are grouped by complete-linkage
#' clustering. A retained group therefore has a maximum pairwise separation
#' no greater than `merge_radius`. Groups containing fewer than
#' `min_cluster_size` input points are discarded.
#'
#' For each retained group, the facet candidate is selected from the original
#' input points rather than from an interpolated centroid. By default, the
#' original point nearest the converged weighted mode is selected;
#' alternatively, the point with the greatest local-height value can be used.
#'
#' @param df A data frame containing thresholded local-height points.
#' @param coord_cols Character vector of length three giving the coordinate
#'   column names. Default: `c("x", "y", "z")`.
#' @param height_col Character. Name of the local-height column used to weight
#'   the condensation. Default: `"height_value"`.
#' @param neighbour_radius Numeric. Radius of the spherical neighbourhood
#'   within which points influence each other during condensation. Must use
#'   the same spatial units as the coordinate columns.
#' @param merge_radius Numeric. Maximum permitted pairwise Euclidean
#'   separation within a final candidate group. Final groups are obtained by
#'   complete-linkage clustering of the converged positions. Default:
#'   `neighbour_radius * 0.6`.
#' @param weight_exponent Numeric, greater than or equal to zero. Exponent
#'   applied to the local-height weights. Non-negative height values are used
#'   directly; values are shifted only when negative heights are present. Larger
#'   values give greater influence to points with high local surface heights.
#'   Default: `2`.
#' @param max_iterations Integer. Maximum number of condensation iterations.
#'   Default: `8`.
#' @param step_size Numeric in the interval `(0, 1]`. Fraction of the
#'   displacement towards the local weighted centre applied during each
#'   iteration. Values below `1` damp the movement. Default: `0.7`.
#' @param tolerance Numeric. Stop iterating when the maximum displacement of
#'   any point during an iteration is less than or equal to this value.
#'   Default: `neighbour_radius * 1e-3`.
#' @param min_cluster_size Integer. Minimum number of original input points
#'   required for a condensed group to be retained as a facet candidate.
#'   Default: `3`.
#' @param select_point Character. Method used to choose the returned candidate
#'   from the original points belonging to each group. `"nearest_mode"`
#'   selects the original point nearest the converged weighted mode;
#'   `"max_height"` selects the point with the greatest local-height value.
#'   Default: `"nearest_mode"`.
#' @param cores Integer. Number of processor cores used for the independent
#'   point-condensation updates. Final grouping is sequential. Default: `1`.
#' @param return_details Logical. If `FALSE`, return only the facet-candidate
#'   data frame. If `TRUE`, additionally return point membership, converged
#'   coordinates, and algorithm parameters. Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress information. Default:
#'   `FALSE`.
#'
#' @return If `return_details = FALSE`, a data frame with one row per retained
#'   facet candidate. Candidate coordinates (`x`, `y`, `z`) correspond to an
#'   original input point. Additional columns report cluster size, selected and
#'   summary height values, converged mode and original centroid coordinates,
#'   the selected source index, and algorithm parameters. If
#'   `return_details = TRUE`, a list with `candidates`, `membership`, and
#'   `parameters`.
#'
#' @export
#' @examples
#' points <- data.frame(
#'   x = c(-0.2, 0, 0.2, 0, 0, 9.8, 10, 10.2, 10, 10),
#'   y = c(0, 0, 0, -0.2, 0.2, 0, 0, 0, -0.2, 0.2),
#'   z = 0,
#'   height_value = c(1, 4, 1, 2, 2, 1, 5, 1, 2, 2)
#' )
#' candidates <- find_facet_candidates_condensed(
#'   points,
#'   neighbour_radius = 0.5,
#'   merge_radius = 0.3,
#'   cores = 1,
#'   verbose = FALSE
#' )
#' candidates
#'
find_facet_candidates_condensed <- function(df,
                                            coord_cols = c("x", "y", "z"),
                                            height_col = "height_value",
                                            neighbour_radius,
                                            merge_radius = neighbour_radius * 0.6,
                                            weight_exponent = 2,
                                            max_iterations = 8,
                                            step_size = 0.7,
                                            tolerance = neighbour_radius * 1e-3,
                                            min_cluster_size = 3,
                                            select_point = c("nearest_mode", "max_height"),
                                            cores = 1,
                                            return_details = FALSE,
                                            verbose = FALSE) {

  select_point <- match.arg(select_point)

  if (!is.data.frame(df)) {
    stop("df must be a data frame.", call. = FALSE)
  }

  if (length(coord_cols) != 3) {
    stop("coord_cols must contain exactly three column names.", call. = FALSE)
  }

  missing_cols <- setdiff(c(coord_cols, height_col), names(df))
  if (length(missing_cols) > 0) {
    stop(
      "df is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (missing(neighbour_radius) || !is.numeric(neighbour_radius) ||
      length(neighbour_radius) != 1 || !is.finite(neighbour_radius) ||
      neighbour_radius <= 0) {
    stop("neighbour_radius must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(merge_radius) || length(merge_radius) != 1 ||
      !is.finite(merge_radius) || merge_radius <= 0) {
    stop("merge_radius must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(weight_exponent) || length(weight_exponent) != 1 ||
      !is.finite(weight_exponent) || weight_exponent < 0) {
    stop("weight_exponent must be a single finite number >= 0.", call. = FALSE)
  }

  if (!is.numeric(step_size) || length(step_size) != 1 ||
      !is.finite(step_size) || step_size <= 0 || step_size > 1) {
    stop("step_size must be in the interval (0, 1].", call. = FALSE)
  }

  if (!is.numeric(max_iterations) || length(max_iterations) != 1 ||
      !is.finite(max_iterations) || max_iterations < 1 ||
      max_iterations != as.integer(max_iterations)) {
    stop("max_iterations must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(tolerance) || length(tolerance) != 1 ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be a finite number greater than or equal to zero.", call. = FALSE)
  }

  if (!is.numeric(min_cluster_size) || length(min_cluster_size) != 1 ||
      !is.finite(min_cluster_size) || min_cluster_size < 1 ||
      min_cluster_size != as.integer(min_cluster_size)) {
    stop("min_cluster_size must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(cores) || length(cores) != 1 ||
      !is.finite(cores) || cores < 1) {
    stop("cores must be a single positive finite number.", call. = FALSE)
  }
  cores <- as.integer(cores)

  empty_candidates <- data.frame(
    facet_candidate_id = integer(),
    x = numeric(),
    y = numeric(),
    z = numeric(),
    cluster_id = integer(),
    point_count = integer(),
    height_value = numeric(),
    height_value_max = numeric(),
    height_value_mean = numeric(),
    mode_x = numeric(),
    mode_y = numeric(),
    mode_z = numeric(),
    centroid_x = numeric(),
    centroid_y = numeric(),
    centroid_z = numeric(),
    selected_source_index = integer(),
    neighbour_radius = numeric(),
    merge_radius = numeric(),
    weight_exponent = numeric(),
    iterations_used = integer(),
    max_shift_final = numeric(),
    cores = integer()
  )

  if (nrow(df) == 0) {
    if (!return_details) {
      return(empty_candidates)
    }
    return(list(
      candidates = empty_candidates,
      membership = data.frame(),
      parameters = list(
        neighbour_radius = neighbour_radius,
        merge_radius = merge_radius,
        weight_exponent = weight_exponent,
        max_iterations = max_iterations,
        step_size = step_size,
        tolerance = tolerance,
        min_cluster_size = min_cluster_size,
        select_point = select_point,
        cores = cores,
        iterations_run = 0L,
        converged = TRUE
      )
    ))
  }

  if (!"source_index" %in% names(df)) {
    df$source_index <- seq_len(nrow(df))
  }

  coords_original <- as.matrix(df[, coord_cols])
  storage.mode(coords_original) <- "double"

  height <- as.numeric(df[[height_col]])

  finite_rows <- stats::complete.cases(coords_original) & is.finite(height)
  if (!any(finite_rows)) {
    stop("No complete finite coordinate/height rows found.", call. = FALSE)
  }

  if (!all(finite_rows)) {
    if (verbose) {
      message("Removing ", sum(!finite_rows), " rows with non-finite coordinates or height values.")
    }
    df <- df[finite_rows, , drop = FALSE]
    coords_original <- coords_original[finite_rows, , drop = FALSE]
    height <- height[finite_rows]
  }

  n <- nrow(df)

  # Use non-negative contrast/height values directly so their absolute
  # 0--1 scale is preserved. Shift only when negative values are present.
  height_min <- min(height, na.rm = TRUE)
  if (height_min < 0) {
    weights <- height - height_min + .Machine$double.eps
  } else {
    weights <- height
  }

  if (weight_exponent != 1) {
    weights <- weights ^ weight_exponent
  }

  if (!any(is.finite(weights)) || sum(weights, na.rm = TRUE) <= 0) {
    weights <- rep(1, n)
  }

  make_grid <- function(coords, cell_size) {
    cell <- floor(coords / cell_size)
    key <- paste(cell[, 1], cell[, 2], cell[, 3], sep = "_")
    grid <- split(seq_len(nrow(coords)), key)
    list(cell = cell, grid = grid)
  }

  neighbour_offsets <- as.matrix(
    expand.grid(
      dx = -1:1,
      dy = -1:1,
      dz = -1:1
    )
  )

  get_neighbours <- function(i, coords, grid_obj, radius_sq) {
    curr_cell <- grid_obj$cell[i, ]

    neighbour_cells <- cbind(
      curr_cell[1] + neighbour_offsets[, 1],
      curr_cell[2] + neighbour_offsets[, 2],
      curr_cell[3] + neighbour_offsets[, 3]
    )

    neighbour_keys <- paste(
      neighbour_cells[, 1],
      neighbour_cells[, 2],
      neighbour_cells[, 3],
      sep = "_"
    )

    candidates <- unlist(grid_obj$grid[neighbour_keys], use.names = FALSE)

    if (length(candidates) == 0) {
      return(i)
    }

    diffs <- sweep(coords[candidates, , drop = FALSE], 2, coords[i, ], "-")
    d2 <- rowSums(diffs * diffs)

    candidates[d2 <= radius_sq]
  }

  coords_current <- coords_original
  neighbour_radius_sq <- neighbour_radius ^ 2
  last_max_shift <- NA_real_
  iterations_used <- 0L

  use_parallel <- cores > 1 && n > 1
  if (use_parallel) {
    cores <- min(cores, n)
    doParallel::registerDoParallel(cores)
    on.exit(doParallel::stopImplicitCluster(), add = TRUE)

    if (verbose) {
      message("Using ", cores, " cores for condensation updates.")
    }
  }

  update_chunk <- function(chunk, coords_current, grid_obj) {
    chunk_out <- coords_current[chunk, , drop = FALSE]

    for (j in seq_along(chunk)) {
      i <- chunk[j]

      neighbours <- get_neighbours(
        i = i,
        coords = coords_current,
        grid_obj = grid_obj,
        radius_sq = neighbour_radius_sq
      )

      w <- weights[neighbours]
      if (!any(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
        next
      }

      local_centre <- colSums(coords_current[neighbours, , drop = FALSE] * w) / sum(w)
      chunk_out[j, ] <- coords_current[i, ] + step_size * (local_centre - coords_current[i, ])
    }

    chunk_out
  }

  for (iter in seq_len(max_iterations)) {
    if (verbose) {
      message("Condensation iteration ", iter, " / ", max_iterations)
    }

    grid_obj <- make_grid(coords_current, neighbour_radius)

    if (use_parallel) {
      # Chunking keeps foreach overhead much lower than launching one task per point.
      n_chunks <- min(n, cores * 4)
      chunk_id <- cut(seq_len(n), breaks = n_chunks, labels = FALSE)
      chunks <- split(seq_len(n), chunk_id)

      chunk <- NULL  # foreach iteration variable; explicit binding for R CMD check
      coords_next <- foreach::foreach(
        chunk = chunks,
        .combine = rbind,
        .inorder = TRUE,
        .export = c("get_neighbours")
      ) %dopar% {
        update_chunk(chunk, coords_current, grid_obj)
      }

    } else {
      coords_next <- coords_current
      coords_next[] <- update_chunk(seq_len(n), coords_current, grid_obj)
    }

    shift <- sqrt(rowSums((coords_next - coords_current) ^ 2))
    last_max_shift <- max(shift, na.rm = TRUE)
    iterations_used <- iter

    coords_current <- coords_next

    if (is.finite(last_max_shift) && last_max_shift <= tolerance) {
      if (verbose) {
        message("Stopping early; max shift = ", signif(last_max_shift, 4))
      }
      break
    }
  }

  if (n == 1L) {
    cluster_id <- 1L
  } else {
    # Complete linkage prevents single-linkage chaining: every pair of
    # converged points within a final group must be no farther apart than
    # merge_radius.
    merge_tree <- stats::hclust(stats::dist(coords_current), method = "complete")
    cluster_id <- stats::cutree(merge_tree, h = merge_radius)
  }

  candidates <- list()
  kept_cluster <- 0L

  for (cl in sort(unique(cluster_id))) {
    idx <- which(cluster_id == cl)

    if (length(idx) < min_cluster_size) {
      next
    }

    kept_cluster <- kept_cluster + 1L

    w <- weights[idx]
    if (!any(is.finite(w)) || sum(w, na.rm = TRUE) <= 0) {
      w <- rep(1, length(idx))
    }

    mode <- colSums(coords_current[idx, , drop = FALSE] * w) / sum(w)
    centroid <- colMeans(coords_original[idx, , drop = FALSE])

    if (select_point == "max_height") {
      selected <- idx[which.max(height[idx])]
    } else {
      d2_to_mode <- rowSums(
        sweep(coords_original[idx, , drop = FALSE], 2, mode, "-") ^ 2
      )

      min_d2 <- min(d2_to_mode, na.rm = TRUE)
      tied <- idx[which(d2_to_mode == min_d2)]

      # If several points are equally close to the mode, prefer the highest one.
      selected <- tied[which.max(height[tied])]
    }

    candidates[[length(candidates) + 1L]] <- data.frame(
      facet_candidate_id = kept_cluster,
      x = coords_original[selected, 1],
      y = coords_original[selected, 2],
      z = coords_original[selected, 3],
      cluster_id = cl,
      point_count = length(idx),
      height_value = height[selected],
      height_value_max = max(height[idx], na.rm = TRUE),
      height_value_mean = mean(height[idx], na.rm = TRUE),
      mode_x = mode[1],
      mode_y = mode[2],
      mode_z = mode[3],
      centroid_x = centroid[1],
      centroid_y = centroid[2],
      centroid_z = centroid[3],
      selected_source_index = df$source_index[selected],
      neighbour_radius = neighbour_radius,
      merge_radius = merge_radius,
      weight_exponent = weight_exponent,
      iterations_used = iterations_used,
      max_shift_final = last_max_shift,
      cores = cores
    )
  }

  if (length(candidates) == 0) {
    candidates_df <- empty_candidates
  } else {
    candidates_df <- do.call(rbind, candidates)
    candidates_df <- candidates_df[order(-candidates_df$height_value_max), , drop = FALSE]
    candidates_df$facet_candidate_id <- seq_len(nrow(candidates_df))
    rownames(candidates_df) <- NULL
  }

  if (!return_details) {
    return(candidates_df)
  }

  membership <- data.frame(
    source_index = df$source_index,
    cluster_id = cluster_id,
    x = coords_original[, 1],
    y = coords_original[, 2],
    z = coords_original[, 3],
    height_value = height,
    converged_x = coords_current[, 1],
    converged_y = coords_current[, 2],
    converged_z = coords_current[, 3]
  )

  list(
    candidates = candidates_df,
    membership = membership,
    parameters = list(
      coord_cols = coord_cols,
      height_col = height_col,
      neighbour_radius = neighbour_radius,
      merge_radius = merge_radius,
      weight_exponent = weight_exponent,
      max_iterations = max_iterations,
      step_size = step_size,
      tolerance = tolerance,
      min_cluster_size = min_cluster_size,
      select_point = select_point,
      iterations_used = iterations_used,
      max_shift_final = last_max_shift,
      cores = cores,
      input_point_count = n,
      candidate_count = nrow(candidates_df)
    )
  )
}