


#' Condense thresholded local-height points into facet candidate modes
#'
#' This function performs a local mode-condensation step on a thresholded
#' local-height point cloud. Each point repeatedly moves toward the
#' height-weighted centre of neighbouring points within a local search radius.
#' After convergence, points that ended close together are merged into one
#' candidate group. For each group, the returned facet candidate is the
#' original input point closest to the converged group mode, not an artificial
#' centroid.
#'
#' This is intended as a pre-clustering or alternative candidate-detection step
#' after local-height thresholding. It avoids a full pairwise distance matrix
#' and is therefore more suitable for thresholded point clouds containing many
#' thousands of points.
#'
#' @param df A data frame containing thresholded local-height points.
#' @param coord_cols Character vector of length 3 giving the coordinate columns.
#'   Defaults to `c("x", "y", "z")`.
#' @param height_col Name of the height/value column used for weighting. For
#'   03B output this is usually `"height_value"`.
#' @param neighbour_radius Radius within which points influence each other
#'   during convergence. Must be in the same units as `x`, `y`, and `z`.
#' @param merge_radius Radius used to merge converged points into final candidate
#'   groups. Defaults to `neighbour_radius / 2`.
#' @param weight_exponent Exponent applied to positive height weights. Higher
#'   values make high local-height points pull more strongly. Defaults to `1`.
#' @param max_iterations Maximum number of condensation iterations. Defaults to
#'   `10`.
#' @param step_size Fraction of the shift toward the local weighted centre used
#'   per iteration. Use values below `1` for more damped convergence. Defaults
#'   to `1`.
#' @param tolerance Stop early when the maximum point displacement in an
#'   iteration falls below this value. Defaults to `neighbour_radius * 1e-3`.
#' @param min_cluster_size Minimum number of original points required for a
#'   converged group to be retained as a facet candidate. Defaults to `1`.
#' @param select_point How to choose the actual returned candidate point within
#'   each converged group. `"nearest_mode"` chooses the original point closest
#'   to the converged group mode. `"max_height"` chooses the highest point in
#'   the group. Defaults to `"nearest_mode"`.
#' @param cores Number of CPU cores used for the independent per-point
#'   condensation updates. The merging step remains sequential. Defaults to `1`.
#' @param return_details Logical. If `FALSE`, return only the facet-candidate
#'   data frame. If `TRUE`, return a list with candidates, membership, converged
#'   coordinates, and parameter metadata.
#' @param verbose Logical. If `TRUE`, print progress messages.
#'
#' @return If `return_details = FALSE`, a data frame with one row per facet
#'   candidate. If `return_details = TRUE`, a list containing candidates,
#'   membership, and parameter metadata.
#'
#' @export
find_facet_candidates_condensed <- function(df,
                                            coord_cols = c("x", "y", "z"),
                                            height_col = "height_value",
                                            neighbour_radius,
                                            merge_radius = neighbour_radius / 2,
                                            weight_exponent = 1,
                                            max_iterations = 10,
                                            step_size = 1,
                                            tolerance = neighbour_radius * 1e-3,
                                            min_cluster_size = 1,
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
    return(empty_candidates)
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

  # Make weights positive while preserving relative height differences.
  height_min <- min(height, na.rm = TRUE)
  weights <- height - height_min + .Machine$double.eps

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

  merge_radius_sq <- merge_radius ^ 2
  grid_final <- make_grid(coords_current, merge_radius)

  cluster_id <- rep(NA_integer_, n)
  curr_cluster <- 0L

  for (i in seq_len(n)) {
    if (!is.na(cluster_id[i])) {
      next
    }

    curr_cluster <- curr_cluster + 1L
    queue <- i
    cluster_id[i] <- curr_cluster

    while (length(queue) > 0) {
      q <- queue[1]
      queue <- queue[-1]

      neighbours <- get_neighbours(
        i = q,
        coords = coords_current,
        grid_obj = grid_final,
        radius_sq = merge_radius_sq
      )

      unassigned <- neighbours[is.na(cluster_id[neighbours])]

      if (length(unassigned) > 0) {
        cluster_id[unassigned] <- curr_cluster
        queue <- c(queue, unassigned)
      }
    }
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