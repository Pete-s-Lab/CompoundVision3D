#' Align a 3D Point Cloud to Anatomical Axes
#'
#' Translates and rotates a three-dimensional point cloud into a standardized
#' landmark-based coordinate system defined by anterior, posterior, left, and
#' right reference landmarks.
#'
#' The anterior landmark is translated to the origin. The positive x-axis is
#' oriented from the right landmark towards the left landmark, and the positive
#' y-axis from the anterior landmark towards the posterior landmark. Because
#' these landmark directions need not be exactly perpendicular, one direction
#' is preserved while the other is orthogonalised according to `priority`. The
#' z-axis is then obtained by a cross product to produce a right-handed
#' coordinate system.
#'
#' A second three-column vector field, normally containing facet surface
#' normals, is rotated by the same transformation. Missing vector triplets are
#' retained as missing values.
#'
#' @param df A data frame or tibble containing point coordinates, landmark
#'   identifiers, and optionally a three-dimensional vector associated with
#'   each point.
#' @param coord_x Character. Name of the x-coordinate column. Default: `"x"`.
#' @param coord_y Character. Name of the y-coordinate column. Default: `"y"`.
#' @param coord_z Character. Name of the z-coordinate column. Default: `"z"`.
#' @param vector_x Character. Name of the x-component of the vectors to rotate.
#'   Default: `"norm.x"`.
#' @param vector_y Character. Name of the y-component of the vectors to rotate.
#'   Default: `"norm.y"`.
#' @param vector_z Character. Name of the z-component of the vectors to rotate.
#'   Default: `"norm.z"`.
#' @param landmark_col Character. Name of the column containing landmark
#'   identifiers. Default: `"ID"`.
#' @param landmark_names Named list defining the values in `landmark_col`
#'   corresponding to the anterior, posterior, left, and right landmarks.
#'   Default: `list(anterior = "anterior", posterior = "posterior",
#'   left = "left", right = "right")`.
#' @param priority Character. Which landmark direction is preserved exactly
#'   during orthogonalisation. `"RL"` preserves the right-to-left axis and
#'   orthogonalises the anterior-to-posterior axis; `"AP"` preserves the
#'   anterior-to-posterior axis and orthogonalises the right-to-left axis.
#'   Default: `"RL"`.
#'
#' @return The input data with the coordinate and vector columns transformed
#'   into the global coordinate system. The returned object additionally has
#'   the attributes `rotation_matrix_rows_xyz`, `translation_applied`,
#'   `priority`, `coordinate_columns`, and `vector_columns`.
#'
#' @export
#' @examples
#' landmarks <- data.frame(
#'   ID = c("anterior", "posterior", "left", "right", "facet_1"),
#'   x = c(0, 0, 1, -1, 0),
#'   y = c(0, 2, 1, 1, 1),
#'   z = c(0, 0, 0, 0, 1),
#'   norm.x = c(NA, NA, NA, NA, 0),
#'   norm.y = c(NA, NA, NA, NA, 0),
#'   norm.z = c(NA, NA, NA, NA, 1)
#' )
#' aligned <- align_pointcloud(landmarks)
#' aligned[, c("ID", "x", "y", "z")]
#'
align_pointcloud <- function(df,
                             coord_x = "x",
                             coord_y = "y",
                             coord_z = "z",
                             vector_x = "norm.x",
                             vector_y = "norm.y",
                             vector_z = "norm.z",
                             landmark_col = "ID",
                             landmark_names = list(
                               anterior = "anterior",
                               posterior = "posterior",
                               left = "left",
                               right = "right"
                             ),
                             priority = "RL") {
  cross3 <- function(a, b) {
    c(a[2] * b[3] - a[3] * b[2],
      a[3] * b[1] - a[1] * b[3],
      a[1] * b[2] - a[2] * b[1])
  }
  nrm <- function(v) {
    s <- sqrt(sum(v^2))
    if (!is.finite(s) || s == 0) stop("Zero-length vector during normalization.")
    v / s
  }
  get_idx <- function(name, lab) {
    idx <- which(lab == name)
    if (length(idx) != 1) {
      stop(sprintf("Expected exactly one '%s' landmark, found %d.", name, length(idx)))
    }
    idx
  }

  coord_cols <- c(coord_x, coord_y, coord_z)
  vector_cols <- c(vector_x, vector_y, vector_z)
  required_cols <- c(coord_cols, landmark_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
  }

  vector_present <- vector_cols %in% names(df)
  if (any(vector_present) && !all(vector_present)) {
    stop("Either all three vector columns or none of them must be present.", call. = FALSE)
  }
  has_vectors <- all(vector_present)
  if (!priority %in% c("RL", "AP")) stop("priority must be 'RL' or 'AP'.")
  if (!all(c("anterior", "posterior", "left", "right") %in% names(landmark_names))) {
    stop("landmark_names must contain named entries: anterior, posterior, left, and right.")
  }

  lab <- df[[landmark_col]]
  ref <- as.matrix(df[, coord_cols])
  storage.mode(ref) <- "double"
  if (any(!is.finite(ref))) {
    stop("Coordinate columns must contain finite values for all rows.")
  }

  iA <- get_idx(landmark_names$anterior, lab)
  iP <- get_idx(landmark_names$posterior, lab)
  iL <- get_idx(landmark_names$left, lab)
  iR <- get_idx(landmark_names$right, lab)

  A <- ref[iA, ]
  P <- ref[iP, ]
  L <- ref[iL, ]
  R <- ref[iR, ]

  Tvec <- as.numeric(A)
  ref_c <- sweep(ref, 2, Tvec, "-")
  v_AP <- as.numeric(P - A)
  v_RL <- as.numeric(L - R)

  if (priority == "AP") {
    y_hat <- nrm(v_AP)
    v_RL_orthY <- v_RL - sum(v_RL * y_hat) * y_hat
    if (sqrt(sum(v_RL_orthY^2)) == 0) {
      stop("Right-Left is collinear with Anterior-Posterior in AP-priority mode.")
    }
    x_hat <- nrm(v_RL_orthY)
    z_hat <- nrm(cross3(x_hat, y_hat))
  } else {
    x_hat <- nrm(v_RL)
    v_AP_orthX <- v_AP - sum(v_AP * x_hat) * x_hat
    if (sqrt(sum(v_AP_orthX^2)) == 0) {
      stop("Anterior-Posterior is collinear with Right-Left in RL-priority mode.")
    }
    y_hat <- nrm(v_AP_orthX)
    z_hat <- nrm(cross3(x_hat, y_hat))
  }

  M <- rbind(x_hat, y_hat, z_hat)

  R_rot <- as.numeric(M %*% (R - A))
  L_rot <- as.numeric(M %*% (L - A))
  if (!(R_rot[1] < L_rot[1])) {
    M[1, ] <- -M[1, ]
    M[3, ] <- -M[3, ]
  }

  ref_rot <- t(M %*% t(ref_c))

  out <- df
  out[[coord_x]] <- ref_rot[, 1]
  out[[coord_y]] <- ref_rot[, 2]
  out[[coord_z]] <- ref_rot[, 3]
  if (has_vectors) {
    vectors <- as.matrix(df[, vector_cols, drop = FALSE])
    storage.mode(vectors) <- "double"
    vectors_rot <- matrix(NA_real_, nrow = nrow(vectors), ncol = 3)
    valid_vectors <- stats::complete.cases(vectors)
    if (any(valid_vectors)) {
      vectors_rot[valid_vectors, ] <- t(M %*% t(vectors[valid_vectors, , drop = FALSE]))
    }
    out[[vector_x]] <- vectors_rot[, 1]
    out[[vector_y]] <- vectors_rot[, 2]
    out[[vector_z]] <- vectors_rot[, 3]
  }

  attr(out, "rotation_matrix_rows_xyz") <- M
  attr(out, "translation_applied") <- Tvec
  attr(out, "priority") <- priority
  attr(out, "coordinate_columns") <- coord_cols
  attr(out, "vector_columns") <- if (has_vectors) vector_cols else character(0)

  out
}


#' Calculate the Normal Vector of a 3D Triangle
#'
#' Calculates a vector perpendicular to the plane defined by three points in
#' three-dimensional space using the cross product of the vectors from `A` to
#' `B` and from `A` to `C`.
#'
#' @param A Numeric vector of length three containing the x, y, and z
#'   coordinates of the first point.
#' @param B Numeric vector of length three containing the x, y, and z
#'   coordinates of the second point.
#' @param C Numeric vector of length three containing the x, y, and z
#'   coordinates of the third point.
#' @param normalize Logical. If `TRUE`, return a unit normal vector. If
#'   `FALSE`, return the unscaled cross product. Default: `TRUE`.
#'
#' @return A numeric vector of length three containing the x, y, and z
#'   components of the triangle normal. If the three points are collinear, a
#'   zero vector is returned.
#'
#' @keywords internal
calculate_normal <- function(A, B, C, normalize = TRUE) {
  A <- as.numeric(A)
  B <- as.numeric(B)
  C <- as.numeric(C)

  if (length(A) != 3 || length(B) != 3 || length(C) != 3 ||
      any(!is.finite(A)) || any(!is.finite(B)) || any(!is.finite(C))) {
    stop("A, B, and C must each contain three finite numeric values.", call. = FALSE)
  }
  if (length(normalize) != 1 || is.na(normalize) || !is.logical(normalize)) {
    stop("normalize must be TRUE or FALSE.", call. = FALSE)
  }

  AB <- B - A
  AC <- C - A
  normal <- c(
    AB[2] * AC[3] - AB[3] * AC[2],
    AB[3] * AC[1] - AB[1] * AC[3],
    AB[1] * AC[2] - AB[2] * AC[1]
  )

  if (normalize) {
    magnitude <- sqrt(sum(normal^2))
    if (magnitude > 0) normal <- normal / magnitude
  }

  normal
}
