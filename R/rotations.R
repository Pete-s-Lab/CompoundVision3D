#' Align 3D Point Cloud to Global Axes
#'
#' Rotates 3D point cloud according to two defined vectors so that these vectors
#' are aligned to the global coordinate system X- and Y-axes.
#'
#' @param df A tibble containing coordinates and IDs.
#' @param ref_x The column containing x-coordinates to be rotated. Default: `x`.
#' @param ref_y The column containing y-coordinates to be rotated. Default: `y`.
#' @param ref_z The column containing z-coordinates to be rotated. Default: `z`.
#' @param data_x  The column containing the reference x-coordinates. Default: `x`.
#' @param data_y  The column containing the reference y-coordinates. Default: `y`.
#' @param data_z  The column containing the reference z-coordinates. Default: `z`.
#' @param landmark_col The column that has to have cells following the 
#' vector defined in `names`. Default: `ID`.
#' @param names A character string defining the global axis to align to.
#' . Default: `anterior, posterior, left, right`. 
#' @param priority A character string defining the global axis that has priority
#' when aligning. Can read `RL` for the right-left vector or `AP` for the 
#' antrior-posterior vector. Default: `RL`.
#' @return Returns a tibble with the aligned coordinates in columns defined by
#' `x_col, y_col, z_col`.
#'
#' @export
#' @examples
#' xxx: add example and change above descsriptionand parameters
#'
align_pointcloud <- function(df,
                             ref_x = "x", 
                             ref_y = "y", 
                             ref_z = "z",           # reference coords (complete for ALL rows)
                             data_x = "norm.x", 
                             data_y = "norm.y", 
                             data_z = "norm.z",     # second coord set (may have NAs anywhere)
                             landmark_col = "ID",
                             names = list(anterior = "anterior",
                                          posterior = "posterior",
                                          left = "left",
                                          right = "right"),
                             priority = "RL"){   # "RL" (right–left first) or "AP" (anterior–posterior first)
  # # testing
  # df = curr_point_cloud
  # ref_x = "x"
  # ref_y = "y"
  # ref_z = "z"
  # data_x = "norm.x"
  # data_y = "norm.y"
  # data_z = "norm.z"
  # landmark_col = "ID"
  # names = list(anterior = "anterior",
  #              posterior = "posterior",
  #              left = "left",
  #              right = "right")
  # priority = "AP"
  
  # helpers
  cross3 <- function(a, b) c(a[2]*b[3]-a[3]*b[2], a[3]*b[1]-a[1]*b[3], a[1]*b[2]-a[2]*b[1])
  nrm <- function(v) { s <- sqrt(sum(v^2)); if (s == 0) stop("Zero-length vector during normalization."); v / s }
  get_idx <- function(name, lab) {
    idx <- which(lab == name)
    if (length(idx) != 1) stop(sprintf("Expected exactly one '%s' landmark, found %d.", name, length(idx)))
    idx
  }
  
  # turn normal vectors into vector end points
  df[, which(colnames(df) == data_x)] <- df[, which(colnames(df) == ref_x)] + 
    df[, which(colnames(df) == data_x)]
  df[, which(colnames(df) == data_y)] <- df[, which(colnames(df) == ref_y)] + 
    df[, which(colnames(df) == data_y)]
  df[, which(colnames(df) == data_z)] <- df[, which(colnames(df) == ref_z)] + 
    df[, which(colnames(df) == data_z)]
  
  lab <- df[[landmark_col]]
  ref <- as.matrix(df[, c(ref_x, ref_y, ref_z)])
  dat <- as.matrix(df[, c(data_x, data_y, data_z)])
  
  if (anyNA(ref))
    stop("Reference coordinates contain NA values; they must be complete for all points.")
  
  iA <- get_idx(names$anterior, lab)
  iP <- get_idx(names$posterior, lab)
  iL <- get_idx(names$left,     lab)
  iR <- get_idx(names$right,    lab)
  
  A <- ref[iA, , drop = FALSE]
  P <- ref[iP, , drop = FALSE]
  L <- ref[iL, , drop = FALSE]
  R <- ref[iR, , drop = FALSE]
  
  # build transform from reference
  Tvec <- as.numeric(A)                 # translate so anterior is origin
  ref_c <- sweep(ref, 2, Tvec, "-")
  
  v_AP <- as.numeric(P - A)             # anterior -> posterior
  v_RL <- as.numeric(L - R)             # right -> left (to satisfy x_right < x_left)
  
  if (priority == "AP") {
    y_hat <- nrm(v_AP)
    v_RL_orthY <- v_RL - sum(v_RL * y_hat) * y_hat
    if (sqrt(sum(v_RL_orthY^2)) == 0) stop("Right–Left is colinear with Anterior–Posterior in AP-priority mode.")
    x_hat <- nrm(v_RL_orthY)
    z_hat <- nrm(cross3(x_hat, y_hat))
  } else if (priority == "RL") {
    x_hat <- nrm(v_RL)
    v_AP_orthX <- v_AP - sum(v_AP * x_hat) * x_hat
    if (sqrt(sum(v_AP_orthX^2)) == 0) stop("Anterior–Posterior is colinear with Right–Left in RL-priority mode.")
    y_hat <- nrm(v_AP_orthX)
    z_hat <- nrm(cross3(x_hat, y_hat))
  } else stop("priority must be 'RL' or 'AP'.")
  
  # rotation matrix; p' = M %*% (p - A)
  M <- rbind(x_hat, y_hat, z_hat)
  
  # enforce orientation on reference landmarks
  R_rot <- as.numeric(M %*% as.numeric(R - A))
  L_rot <- as.numeric(M %*% as.numeric(L - A))
  if (!(R_rot[1] < L_rot[1])) {
    M[1, ] <- -M[1, ]; R_rot[1] <- -R_rot[1]; L_rot[1] <- -L_rot[1]
  }
  if (!(0 < R_rot[2])) {
    M[2, ] <- -M[2, ]; M[3, ] <- -M[3, ]  # keep right-handedness
  }
  
  # apply to BOTH coordinate sets
  # reference (complete): rotate all rows
  ref_rot <- t(M %*% t(ref_c))
  
  # data (may have NA anywhere): only rotate complete triplets
  dat_c <- sweep(dat, 2, Tvec, "-")
  dat_rot <- matrix(NA_real_, nrow = nrow(dat_c), ncol = 3)
  valid <- stats::complete.cases(dat_c)
  if (any(valid)) dat_rot[valid, ] <- t(M %*% t(dat_c[valid, , drop = FALSE]))
  
  # write back
  out <- df
  out[[ref_x]]  <- ref_rot[, 1]; out[[ref_y]] <- ref_rot[, 2]; out[[ref_z]] <- ref_rot[, 3]
  out[[data_x]] <- dat_rot[, 1]; out[[data_y]] <- dat_rot[, 2]; out[[data_z]] <- dat_rot[, 3]
  
  # turn vector end points into normal vectors
  out[, which(colnames(out) == data_x)] <- out[, which(colnames(out) == data_x)] - 
    out[, which(colnames(out) == ref_x)]
  out[, which(colnames(out) == data_y)] <- out[, which(colnames(out) == data_y)] - 
    out[, which(colnames(out) == ref_y)]
  out[, which(colnames(out) == data_z)] <- out[, which(colnames(out) == data_z)] -
    out[, which(colnames(out) == ref_z)]
  
  
  # metadata
  attr(out, "rotation_matrix_rows_xyz") <- M
  attr(out, "translation_applied") <- Tvec
  attr(out, "priority") <- priority
  attr(out, "reference_columns") <- c(ref_x, ref_y, ref_z)
  attr(out, "data_columns")      <- c(data_x, data_y, data_z)
  
  return(out)
}




#' Calculate normal of triangle in 3D
#'
#' xxx: add description
#'
#' @param A A `vector` with the `numeric` x, y, and z coordinges of point 1.
#' @param B A `vector` with the `numeric` x, y, and z coordinges of point 2.
#' @param C A `vector` with the `numeric` x, y, and z coordinges of point 3.
#' @return Returns a `vector` with three `numeric` data points describing the 
#' normal vector of the triangle.
#'
#' @export
#' @examples
#' xxx: add example
#'
calculate_normal <- function(A, B, C, normalize = TRUE) {
  # Ensure the points are numeric vectors of length 3
  if (length(A) != 3 || length(B) != 3 || length(C) != 3) {
    stop("All points must be numeric vectors of length 3.")
  }
  
  # Calculate the vectors AB and AC
  AB <- B - A
  AC <- C - A
  
  # Compute the cross product AB x AC
  normal <- c(
    AB[2] * AC[3] - AB[3] * AC[2],
    AB[3] * AC[1] - AB[1] * AC[3],
    AB[1] * AC[2] - AB[2] * AC[1]
  )
  
  # Calculate the magnitude of the normal vector
  magnitude <- sqrt(sum(normal^2))
  
  # Normalize the normal vector to get a unit normal vector
  if(normalize == TRUE){
    if (magnitude != 0) {
      normal <- normal / magnitude
    }
  }
  
  return(normal)
}


#' Calculate angle between two vectors in 3D
#'
#' xxx: add description
#'
#' @param a A `vector` with three `numeric` data points describing vector 1.
#' @param b A `vector` with three `numeric` data points describing vector b.
#' @return Returns a `numeric` value with the angle in degree (°).
#'
#' @export
#' @examples
#' xxx: add example
#'
angle_between_vectors <- function(a, b) {
  # # testing
  # a = curr_normal
  # b = center_vector
  
  # # plot vectors starting at 0
  # # plot3d(x=NULL)
  # lines3d(x = c(0, curr_normal[1]),
  #         y = c(0, curr_normal[2]),
  #         z = c(0, curr_normal[3]),
  #         col = "cyan")
  # lines3d(x = c(0, center_vector[1]),
  #         y = c(0, center_vector[2]),
  #         z = c(0, center_vector[3]),
  #         col = "red")
  
  # Compute the dot product
  dot_product <- sum(a * b)
  
  # Compute the magnitudes of the vectors
  magnitude_a <- sqrt(sum(a^2))
  magnitude_b <- sqrt(sum(b^2))
  
  # Compute the cosine of the angle between the vectors
  cos_theta <- dot_product / (magnitude_a * magnitude_b)
  angle_deg <- cos_theta*180/pi
  return(angle_deg)
}