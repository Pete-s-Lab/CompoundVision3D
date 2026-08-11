# Internal state used to retain the most recent face-on 2D and 3D views.
# This allows add = TRUE to redraw/extend a 2D comparison plot and to reuse
# the reference orientation and dimensions in 3D. It is intentionally not
# exported.
.cv3d_face_on_state <- new.env(parent = emptyenv())

#' Plot a Compound Eye in a Standardised Face-On View
#'
#' Displays a compound-eye point cloud approximately face-on without modifying
#' the input data. The viewing direction is estimated from the mean surface-
#' normal direction. The eye can be displayed either as a two-dimensional
#' base-R projection or as an interactively rotatable three-dimensional `rgl`
#' scene.
#'
#' Each cloud is centred independently for display. Optional translations are
#' therefore deliberate plotting offsets rather than translations inherited
#' from the input coordinate system.
#'
#' @param df A data frame or tibble containing three-dimensional point
#'   coordinates and, unless `view` or `view_direction` is supplied,
#'   surface-normal vectors.
#' @param coord_x Character. Name of the x-coordinate column. Default: `"x"`.
#' @param coord_y Character. Name of the y-coordinate column. Default: `"y"`.
#' @param coord_z Character. Name of the z-coordinate column. Default: `"z"`.
#' @param normal_x Character. Name of the x-component of the surface normals.
#'   Default: `"norm.x"`.
#' @param normal_y Character. Name of the y-component of the surface normals.
#'   Default: `"norm.y"`.
#' @param normal_z Character. Name of the z-component of the surface normals.
#'   Default: `"norm.z"`.
#' @param projection Character. Plotting mode. `"2D"` produces a face-on
#'   orthographic projection using base R. `"3D"` displays the centred cloud
#'   using `rgl` and applies the equivalent face-on viewpoint. Default:
#'   `c("2D", "3D")`, i.e. `"2D"`.
#' @param view Optional view object previously returned by
#'   `view_eye_face_on()`. Supplying it reuses the same viewing direction,
#'   in-plane orientation, and reference dimensions. Default: `NULL`.
#' @param view_direction Optional numeric vector of length three specifying the
#'   outward viewing direction. If `NULL`, it is estimated from the mean of the
#'   valid unit surface normals. Ignored when `view` is supplied.
#' @param up Numeric vector of length three giving the preferred vertical
#'   direction before automatic in-plane orientation. Default: `c(0, 1, 0)`.
#' @param reverse Logical. If `TRUE`, reverse the inferred viewing direction.
#'   Ignored when `view` is supplied. Default: `FALSE`.
#' @param long_axis_vertical Logical. If `TRUE`, rotate the face-on view by
#'   90 degrees when necessary so that the greater projected coordinate range
#'   of the reference cloud lies along the displayed y-axis. No PCA is used.
#'   Ignored when `view` is supplied. Default: `TRUE`.
#' @param add Logical. If `FALSE`, create a new plot or `rgl` scene. If `TRUE`,
#'   add a new comparison cloud to the most recent CV3D face-on plot of the
#'   same projection type. Default: `FALSE`.
#' @param translate_x Numeric. Horizontal display translation in units of the
#'   face-on width of the first/reference cloud. Thus `translate_x = 1` moves
#'   the current cloud by exactly one reference-cloud width. Default: `0`.
#' @param translate_y Numeric. Vertical display translation in units of the
#'   face-on height of the first/reference cloud. Thus `translate_y = 1` moves
#'   the current cloud by exactly one reference-cloud height. Default: `0`.
#' @param col Point colour or vector of point colours. Default: `"black"`.
#' @param pch Plotting symbol for the two-dimensional projection. The default
#'   `"."` is suitable for dense surface point clouds.
#' @param cex Point-size multiplier for the two-dimensional projection.
#'   Default: `1`.
#' @param rgl_size Point size for three-dimensional `rgl` plotting.
#'   Default: `2`.
#' @param fov Numeric field-of-view angle for the `rgl` view. `0` gives an
#'   orthographic projection. Default: `0`.
#' @param zoom Numeric zoom factor for the `rgl` view. Default: `1`.
#' @param axes Logical. Whether axes should be drawn for a new plot.
#'   Default: `TRUE`.
#' @param xlab Character. x-axis label for a new 2D plot. Default: `"x"`.
#' @param ylab Character. y-axis label for a new 2D plot. Default: `"y"`.
#'
#' @details
#' Each valid surface normal is first normalised to unit length so that vector
#' magnitude does not influence the estimated viewing direction. Their mean
#' direction is then normalised and used as the viewing axis.
#'
#' The preferred `up` vector is projected onto the plane perpendicular to the
#' viewing direction. Together with the viewing direction, this defines the
#' face-on coordinate system. If `long_axis_vertical = TRUE`, only a possible
#' 90-degree in-plane rotation is added so that the greater projected range of
#' the reference eye runs vertically. No principal-component analysis is
#' performed.
#'
#' Before plotting, every supplied cloud is centred on its own arithmetic
#' coordinate centroid. Any original translation between separately supplied
#' clouds is therefore intentionally discarded. `translate_x` and
#' `translate_y` then provide explicit display offsets relative to the width
#' and height of the first/reference cloud. For example, `translate_x = 1.1`
#' places a second similarly sized eye to the right of the first with a gap of
#' approximately 10 percent of the reference width.
#'
#' With `projection = "2D"`, `add = TRUE` redraws the CV3D-managed layers so
#' that the plotting limits expand automatically when newly added clouds fall
#' outside the previous range. Because base graphics are redrawn, unrelated
#' annotations added manually between calls are not retained by this redraw.
#'
#' With `projection = "3D"`, the `rgl` package is required only at run time and
#' can therefore remain in `Suggests`. Added clouds use the same face-on camera
#' and the same translation convention as the 2D projection.
#'
#' CV3D assumes spatial coordinates are supplied in micrometres (µm).
#'
#' @return Invisibly returns a list containing:
#'
#' * `view_direction`: unit viewing direction in the original coordinate
#'   system;
#' * `x_direction`: horizontal screen direction in the original coordinate
#'   system;
#' * `y_direction`: vertical screen direction in the original coordinate
#'   system;
#' * `reference_width`: face-on width of the first/reference cloud;
#' * `reference_height`: face-on height of the first/reference cloud;
#' * `userMatrix`: 4 by 4 matrix used for the equivalent `rgl` view;
#' * `rotated_90`: whether the automatic 90-degree in-plane rotation was
#'   applied;
#' * `normal_resultant_length`: length of the mean unit-normal vector;
#' * `cloud_centre`: arithmetic coordinate centroid of the current cloud;
#' * `translate_x` and `translate_y`: applied dimensionless display offsets;
#' * `translation_3d`: corresponding three-dimensional display translation;
#'   and
#' * `projected_coordinates`: two-dimensional coordinates actually plotted.
#'
#' The input data are not modified.
#'
#' @seealso [align_pointcloud()]
#' @export
#'
#' @examples
#' data(cv3d_example_surface)
#'
#' # Standardised base-R face-on projection; rgl is not required.
#' eye_view <- view_eye_face_on(
#'   cv3d_example_surface,
#'   projection = "2D"
#' )
#'
#' # Add a second centred copy 1.1 reference-eye widths to the right.
#' view_eye_face_on(
#'   cv3d_example_surface,
#'   projection = "2D",
#'   add = TRUE,
#'   translate_x = 1.1,
#'   col = "grey50"
#' )
#'
#' if (interactive() && requireNamespace("rgl", quietly = TRUE)) {
#'   view_eye_face_on(
#'     cv3d_example_surface,
#'     projection = "3D"
#'   )
#'
#'   view_eye_face_on(
#'     cv3d_example_surface,
#'     projection = "3D",
#'     add = TRUE,
#'     translate_x = 1.1,
#'     col = "grey50"
#'   )
#' }

view_eye_face_on <- function(df,
                             coord_x = "x",
                             coord_y = "y",
                             coord_z = "z",
                             normal_x = "norm.x",
                             normal_y = "norm.y",
                             normal_z = "norm.z",
                             projection = c("2D", "3D"),
                             view = NULL,
                             view_direction = NULL,
                             up = c(0, 1, 0),
                             reverse = FALSE,
                             long_axis_vertical = TRUE,
                             add = FALSE,
                             translate_x = 0,
                             translate_y = 0,
                             col = "black",
                             pch = ".",
                             cex = 1,
                             rgl_size = 2,
                             fov = 0,
                             zoom = 1,
                             axes = TRUE,
                             xlab = "x",
                             ylab = "y") {
  
  state <- .cv3d_face_on_state
  
  projection <- match.arg(projection)
  
  if (length(translate_x) != 1 ||
      !is.finite(translate_x) ||
      length(translate_y) != 1 ||
      !is.finite(translate_y)) {
    stop(
      "`translate_x` and `translate_y` must be single finite numbers.",
      call. = FALSE
    )
  }
  
  cross3 <- function(a, b) {
    c(
      a[2] * b[3] - a[3] * b[2],
      a[3] * b[1] - a[1] * b[3],
      a[1] * b[2] - a[2] * b[1]
    )
  }
  
  unit3 <- function(v, label) {
    
    v <- as.numeric(v)
    
    if (length(v) != 3 || any(!is.finite(v))) {
      stop(
        label,
        " must contain three finite numeric values.",
        call. = FALSE
      )
    }
    
    len <- sqrt(sum(v^2))
    
    if (!is.finite(len) ||
        len <= sqrt(.Machine$double.eps)) {
      stop(
        label,
        " must have non-zero length.",
        call. = FALSE
      )
    }
    
    v / len
  }
  
  coord_cols <- c(coord_x, coord_y, coord_z)
  missing_coords <- setdiff(coord_cols, names(df))
  
  if (length(missing_coords) > 0) {
    stop(
      "Missing coordinate column(s): ",
      paste(missing_coords, collapse = ", "),
      call. = FALSE
    )
  }
  
  coords <- as.matrix(
    df[, coord_cols, drop = FALSE]
  )
  
  storage.mode(coords) <- "double"
  
  if (nrow(coords) < 1 ||
      any(!is.finite(coords))) {
    stop(
      "Coordinate columns must contain finite values.",
      call. = FALSE
    )
  }
  
  normal_resultant_length <- NA_real_
  rotated_90 <- FALSE
  
  state_name <- if (projection == "2D") {
    "view_2d"
  } else {
    "view_3d"
  }
  
  # ------------------------------------------------------------
  # Reuse the previous view for added layers unless one is
  # supplied explicitly.
  # ------------------------------------------------------------
  
  if (add && is.null(view)) {
    
    if (!exists(
      state_name,
      envir = state,
      inherits = FALSE
    )) {
      stop(
        paste0(
          "No existing ",
          projection,
          " face-on plot is available. ",
          "Create one with add = FALSE first."
        ),
        call. = FALSE
      )
    }
    
    view <- get(
      state_name,
      envir = state,
      inherits = FALSE
    )$view
  }
  
  # ------------------------------------------------------------
  # Determine or reuse the orientation.
  # ------------------------------------------------------------
  
  if (!is.null(view)) {
    
    required <- c(
      "view_direction",
      "x_direction",
      "y_direction",
      "reference_width",
      "reference_height"
    )
    
    if (!all(required %in% names(view))) {
      stop(
        "`view` is not a valid result from this version of ",
        "view_eye_face_on().",
        call. = FALSE
      )
    }
    
    view_direction <-
      as.numeric(view$view_direction)
    
    x_screen <-
      as.numeric(view$x_direction)
    
    y_screen <-
      as.numeric(view$y_direction)
    
    reference_width <-
      as.numeric(view$reference_width)
    
    reference_height <-
      as.numeric(view$reference_height)
    
    rotated_90 <-
      isTRUE(view$rotated_90)
    
    normal_resultant_length <-
      view$normal_resultant_length
    
  } else {
    
    if (is.null(view_direction)) {
      
      normal_cols <- c(
        normal_x,
        normal_y,
        normal_z
      )
      
      missing_normals <- setdiff(
        normal_cols,
        names(df)
      )
      
      if (length(missing_normals) > 0) {
        stop(
          "Missing surface-normal column(s): ",
          paste(missing_normals, collapse = ", "),
          call. = FALSE
        )
      }
      
      normals <- as.matrix(
        df[, normal_cols, drop = FALSE]
      )
      
      storage.mode(normals) <- "double"
      
      lengths <- sqrt(rowSums(normals^2))
      
      valid <- is.finite(lengths) &
        lengths > sqrt(.Machine$double.eps)
      
      if (!any(valid)) {
        stop(
          "No valid non-zero surface normals are available.",
          call. = FALSE
        )
      }
      
      unit_normals <-
        normals[valid, , drop = FALSE] /
        lengths[valid]
      
      mean_normal <- colMeans(unit_normals)
      
      normal_resultant_length <-
        sqrt(sum(mean_normal^2))
      
      if (!is.finite(normal_resultant_length) ||
          normal_resultant_length <=
          sqrt(.Machine$double.eps)) {
        
        stop(
          paste(
            "The surface normals do not define",
            "a clear mean direction."
          ),
          call. = FALSE
        )
      }
      
      view_direction <-
        unname(
          mean_normal /
            normal_resultant_length
        )
      
    } else {
      
      view_direction <- unit3(
        view_direction,
        "view_direction"
      )
    }
    
    if (reverse) {
      view_direction <- -view_direction
    }
    
    up <- unit3(up, "up")
    
    up_projected <-
      up -
      sum(up * view_direction) *
      view_direction
    
    up_length <-
      sqrt(sum(up_projected^2))
    
    if (up_length <= 1e-8) {
      
      candidate_axes <- diag(3)
      
      fallback <- candidate_axes[
        which.min(
          abs(
            as.numeric(
              candidate_axes %*%
                view_direction
            )
          )
        ),
      ]
      
      up_projected <-
        fallback -
        sum(fallback * view_direction) *
        view_direction
      
      up_length <-
        sqrt(sum(up_projected^2))
    }
    
    y_screen <-
      up_projected / up_length
    
    x_screen <- cross3(
      y_screen,
      view_direction
    )
    
    x_screen <-
      unname(
        x_screen /
          sqrt(sum(x_screen^2))
      )
    
    y_screen <- cross3(
      view_direction,
      x_screen
    )
    
    y_screen <-
      unname(
        y_screen /
          sqrt(sum(y_screen^2))
      )
    
    # Determine orientation from the first/reference cloud.
    cloud_centre <- colMeans(coords)
    
    coords_centred <- sweep(
      coords,
      2,
      cloud_centre,
      "-"
    )
    
    projected_x <- as.numeric(
      coords_centred %*% x_screen
    )
    
    projected_y <- as.numeric(
      coords_centred %*% y_screen
    )
    
    if (long_axis_vertical &&
        diff(range(projected_x)) >
        diff(range(projected_y))) {
      
      old_x_screen <- x_screen
      
      x_screen <- unname(-y_screen)
      y_screen <- unname(old_x_screen)
      
      projected_x <- as.numeric(
        coords_centred %*% x_screen
      )
      
      projected_y <- as.numeric(
        coords_centred %*% y_screen
      )
      
      rotated_90 <- TRUE
    }
    
    reference_width <-
      diff(range(projected_x))
    
    reference_height <-
      diff(range(projected_y))
    
    if (reference_width <= 0 ||
        reference_height <= 0) {
      stop(
        "The reference cloud must span both displayed dimensions.",
        call. = FALSE
      )
    }
  }
  
  # ------------------------------------------------------------
  # Centre THIS cloud independently.
  # Original specimen translation is deliberately discarded.
  # ------------------------------------------------------------
  
  cloud_centre <- colMeans(coords)
  
  coords_centred <- sweep(
    coords,
    2,
    cloud_centre,
    "-"
  )
  
  projected_x <- as.numeric(
    coords_centred %*% x_screen
  )
  
  projected_y <- as.numeric(
    coords_centred %*% y_screen
  )
  
  # ------------------------------------------------------------
  # User-defined display translation.
  #
  # translate_x = 1:
  # one reference-eye width to the right.
  #
  # translate_y = 1:
  # one reference-eye height upwards.
  # ------------------------------------------------------------
  
  shift_x <-
    translate_x * reference_width
  
  shift_y <-
    translate_y * reference_height
  
  projected_x <-
    projected_x + shift_x
  
  projected_y <-
    projected_y + shift_y
  
  # Same translation expressed in original 3D coordinates.
  translation_3d <-
    unname(
      shift_x * x_screen +
        shift_y * y_screen
    )
  
  coords_plot <- sweep(
    coords_centred,
    2,
    translation_3d,
    "+"
  )
  
  user_matrix <- diag(4)
  
  user_matrix[1, 1:3] <-
    x_screen
  
  user_matrix[2, 1:3] <-
    y_screen
  
  user_matrix[3, 1:3] <-
    view_direction
  
  result_view <- list(
    view_direction =
      unname(view_direction),
    x_direction =
      unname(x_screen),
    y_direction =
      unname(y_screen),
    reference_width =
      reference_width,
    reference_height =
      reference_height,
    userMatrix =
      user_matrix,
    rotated_90 =
      rotated_90,
    normal_resultant_length =
      normal_resultant_length
  )
  
  result <- c(
    result_view,
    list(
      cloud_centre =
        cloud_centre,
      translate_x =
        translate_x,
      translate_y =
        translate_y,
      translation_3d =
        translation_3d,
      projected_coordinates =
        cbind(
          x = projected_x,
          y = projected_y
        )
    )
  )
  
  # ============================================================
  # 2D
  # ============================================================
  
  if (projection == "2D") {
    
    layer <- list(
      x = projected_x,
      y = projected_y,
      col = col,
      pch = pch,
      cex = cex
    )
    
    draw_layers <- function(plot_state) {
      
      all_x <- unlist(
        lapply(
          plot_state$layers,
          `[[`,
          "x"
        ),
        use.names = FALSE
      )
      
      all_y <- unlist(
        lapply(
          plot_state$layers,
          `[[`,
          "y"
        ),
        use.names = FALSE
      )
      
      graphics::plot(
        NA_real_,
        NA_real_,
        type = "n",
        xlim = range(all_x),
        ylim = range(all_y),
        asp = 1,
        axes = plot_state$axes,
        xlab = plot_state$xlab,
        ylab = plot_state$ylab
      )
      
      for (current_layer in
           plot_state$layers) {
        
        graphics::points(
          current_layer$x,
          current_layer$y,
          col = current_layer$col,
          pch = current_layer$pch,
          cex = current_layer$cex
        )
      }
    }
    
    if (!add) {
      
      plot_state <- list(
        view = result_view,
        layers = list(layer),
        axes = axes,
        xlab = xlab,
        ylab = ylab
      )
      
      draw_layers(plot_state)
      
      assign(
        state_name,
        plot_state,
        envir = state
      )
      
    } else {
      
      plot_state <- get(
        state_name,
        envir = state,
        inherits = FALSE
      )
      
      plot_state$layers[[length(plot_state$layers) + 1]] <- layer
      
      draw_layers(plot_state)
      
      assign(
        state_name,
        plot_state,
        envir = state
      )
    }
    
    # ============================================================
    # 3D
    # ============================================================
    
  } else {
    
    if (!requireNamespace(
      "rgl",
      quietly = TRUE
    )) {
      stop(
        'Package "rgl" is required for projection = "3D".',
        call. = FALSE
      )
    }
    
    if (!add) {
      
      rgl::open3d()
      
      rgl::plot3d(
        coords_plot[, 1],
        coords_plot[, 2],
        coords_plot[, 3],
        type = "p",
        col = col,
        size = rgl_size,
        aspect = "iso",
        axes = axes,
        xlab = coord_x,
        ylab = coord_y,
        zlab = coord_z
      )
      
      rgl::view3d(
        userMatrix = user_matrix,
        fov = fov,
        zoom = zoom,
        scale = c(1, 1, 1)
      )
      
      assign(
        state_name,
        list(
          view = result_view
        ),
        envir = state
      )
      
    } else {
      
      if (rgl::cur3d() == 0) {
        stop(
          "No active rgl device is available for add = TRUE.",
          call. = FALSE
        )
      }
      
      rgl::points3d(
        coords_plot[, 1],
        coords_plot[, 2],
        coords_plot[, 3],
        col = col,
        size = rgl_size
      )
      
      rgl::view3d(
        userMatrix = user_matrix,
        fov = fov,
        zoom = zoom,
        scale = c(1, 1, 1)
      )
    }
  }
  
  invisible(result)
}
