#' Align 3D Point Cloud to Global Axes
#'
#' Rotates 3D point cloud according to two defined vectors so that these vectors
#' are aligned to the global coordinate system X- and Y-axes.
#'
#' @param df A tibble containing coordinates and IDs.
#' @param x_col The column containing x-coordinates of the points. Default: `x`.
#' @param y_col The column containing y-coordinates of the points. Default: `y`.
#' @param z_col The column containing z-coordinates of the points. Default: `z`.
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
align_pointcloud <- function(
    df,
    x_col = "x", y_col = "y", z_col = "z",
    landmark_col = "ID",
    names = list(anterior = "anterior",
                 posterior = "posterior",
                 left = "left",
                 right = "right"),
    priority = "RL") { # "RL" (right–left first) or "AP" (anterior–posterior first)
  # priority <- match.arg(priority)
  
  # Helpers
  cross3 <- function(a, b) c(
    a[2]*b[3] - a[3]*b[2],
    a[3]*b[1] - a[1]*b[3],
    a[1]*b[2] - a[2]*b[1]
  )
  nrm <- function(v) {
    s <- sqrt(sum(v^2))
    if (s == 0) stop("Zero-length vector during normalization.")
    v / s
  }
  
  # Pull coords
  xyz <- as.matrix(df[, c(x_col, y_col, z_col)])
  if (!all(is.finite(xyz))) stop("Non-finite coordinates found.")
  lab <- df[[landmark_col]]
  
  get_pt <- function(name) {
    idx <- which(lab == name)
    if (length(idx) != 1) stop(sprintf("Expected exactly one '%s' landmark, found %d.", name, length(idx)))
    xyz[idx, , drop = FALSE]
  }
  A <- get_pt(names$anterior)
  P <- get_pt(names$posterior)
  L <- get_pt(names$left)
  R <- get_pt(names$right)
  
  # Translate so anterior is the origin (convenient for the y_anterior < y_right check)
  Tvec <- as.numeric(A)
  xyz_c <- sweep(xyz, 2, Tvec, "-")
  
  # Base vectors
  v_AP <- as.numeric(P - A)  # anterior -> posterior
  v_RL <- as.numeric(L - R)  # right -> left (note: we choose R->L to match the final x_right < x_left)
  # (If you prefer L->R, it also works; we'll enforce orientation later.)
  
  # Build orthonormal axes based on priority
  if (priority == "AP") {
    # 1) Y along AP
    y_hat <- nrm(v_AP)
    # 2) X from RL projected off Y
    v_RL_orthY <- v_RL - sum(v_RL * y_hat) * y_hat
    if (sqrt(sum(v_RL_orthY^2)) == 0)
      stop("Right–Left is colinear with Anterior–Posterior; cannot define X in AP-priority mode.")
    x_hat <- nrm(v_RL_orthY)
    # 3) Right-handed Z
    z_hat <- nrm(cross3(x_hat, y_hat))
  } else { # priority == "RL"
    # 1) X along RL
    x_hat <- nrm(v_RL)
    # 2) Y from AP projected off X
    v_AP_orthX <- v_AP - sum(v_AP * x_hat) * x_hat
    if (sqrt(sum(v_AP_orthX^2)) == 0)
      stop("Anterior–Posterior is colinear with Right–Left; cannot define Y in RL-priority mode.")
    y_hat <- nrm(v_AP_orthX)
    # 3) Right-handed Z
    z_hat <- nrm(cross3(x_hat, y_hat))
  }
  
  # Rotation: rows are new axes (X, Y, Z) expressed in original basis
  # p' = M %*% (p - A)
  M <- rbind(x_hat, y_hat, z_hat)
  
  # Apply rotation
  xyz_r <- t(M %*% t(xyz_c))
  
  # --- Enforce required orientations ---------------------------------------
  
  # 1) Enforce x_right < x_left
  R_rot <- as.numeric(M %*% as.numeric(R - A))
  L_rot <- as.numeric(M %*% as.numeric(L - A))
  if (!(R_rot[1] < L_rot[1])) {
    # Flip X axis
    M[1, ] <- -M[1, ]
    xyz_r[, 1] <- -xyz_r[, 1]
    # update rotated landmarks for subsequent checks
    R_rot[1] <- -R_rot[1]
    L_rot[1] <- -L_rot[1]
  }
  
  # 2) Enforce y_anterior < y_right
  # After translation, anterior is at (0,0,0) so condition is simply 0 < y(R)
  if (!(0 < R_rot[2])) {
    # Flip Y (and Z to keep right-handedness)
    M[2, ] <- -M[2, ]
    M[3, ] <- -M[3, ]
    xyz_r[, 2] <- -xyz_r[, 2]
    xyz_r[, 3] <- -xyz_r[, 3]
  }
  
  # --------------------------------------------------------------------------
  
  out <- df
  out[[x_col]] <- xyz_r[, 1]
  out[[y_col]] <- xyz_r[, 2]
  out[[z_col]] <- xyz_r[, 3]
  
  attr(out, "rotation_matrix_rows_xyz") <- M
  attr(out, "translation_applied") <- Tvec
  attr(out, "priority") <- priority
  out
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