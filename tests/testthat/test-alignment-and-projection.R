test_that("align_pointcloud establishes the CV3D anatomical coordinate system", {
  landmarks <- data.frame(
    ID = c("anterior", "posterior", "left", "right", "facet_1"),
    x = c(2, 2, 3, 1, 2),
    y = c(3, 5, 4, 4, 4),
    z = c(4, 4, 4, 4, 5),
    norm.x = c(NA, NA, NA, NA, 0),
    norm.y = c(NA, NA, NA, NA, 0),
    norm.z = c(NA, NA, NA, NA, 1)
  )

  result <- align_pointcloud(landmarks, priority = "RL")

  anterior <- result[result$ID == "anterior", c("x", "y", "z")]
  posterior <- result[result$ID == "posterior", c("x", "y", "z")]
  left <- result[result$ID == "left", c("x", "y", "z")]
  right <- result[result$ID == "right", c("x", "y", "z")]

  expect_equal(as.numeric(anterior), c(0, 0, 0), tolerance = 1e-12)
  expect_gt(posterior$y, anterior$y)
  expect_gt(left$x, right$x)

  M <- attr(result, "rotation_matrix_rows_xyz")
  expect_equal(unname(M %*% t(M)), diag(3), tolerance = 1e-12)
  expect_equal(det(M), 1, tolerance = 1e-12)
})


test_that("align_pointcloud works without an associated vector field", {
  landmarks <- data.frame(
    ID = c("anterior", "posterior", "left", "right", "facet_1"),
    x = c(2, 2, 3, 1, 2),
    y = c(3, 5, 4, 4, 4),
    z = c(4, 4, 4, 4, 5)
  )

  result <- align_pointcloud(landmarks, priority = "RL")
  expect_false(any(c("norm.x", "norm.y", "norm.z") %in% names(result)))
  expect_equal(attr(result, "vector_columns"), character(0))
})

test_that("align_pointcloud rotates vectors without applying translation", {
  landmarks <- data.frame(
    ID = c("anterior", "posterior", "left", "right", "facet_1"),
    x = c(2, 2, 3, 1, 2),
    y = c(3, 5, 4, 4, 4),
    z = c(4, 4, 4, 4, 5),
    norm.x = c(NA, NA, NA, NA, 1),
    norm.y = c(NA, NA, NA, NA, 0),
    norm.z = c(NA, NA, NA, NA, 0)
  )

  result <- align_pointcloud(landmarks, priority = "RL")
  M <- attr(result, "rotation_matrix_rows_xyz")
  expected <- as.numeric(M %*% c(1, 0, 0))
  got <- as.numeric(result[result$ID == "facet_1", c("norm.x", "norm.y", "norm.z")])
  expect_equal(got, expected, tolerance = 1e-12)
})


test_that("ray_sphere_intersection handles forward, tangent and absent intersections", {
  expect_equal(
    ray_sphere_intersection(c(0, 0, 0), c(1, 0, 0), c(0, 0, 0), 10),
    c(x = 10, y = 0, z = 0),
    tolerance = 1e-12
  )

  expect_equal(
    ray_sphere_intersection(c(-1, 1, 0), c(1, 0, 0), c(0, 0, 0), 1),
    c(x = 0, y = 1, z = 0),
    tolerance = 1e-12
  )

  no_hit <- ray_sphere_intersection(c(2, 0, 0), c(1, 0, 0), c(0, 0, 0), 1)
  expect_true(all(is.na(no_hit)))

  expect_equal(
    ray_sphere_intersection(c(5, 0, 0), c(0, 1, 0), c(5, 0, 0), 2),
    c(x = 5, y = 2, z = 0),
    tolerance = 1e-12
  )
})


test_that("cartesian_to_view_angles follows the CV3D viewing convention", {
  result <- cartesian_to_view_angles(
    x = c(0, 1, 0, -1, 0),
    y = c(-1, 0, 1, 0, 0),
    z = c(0, 0, 0, 0, 1)
  )

  expect_equal(result$azimuth[1:4], c(0, 90, -180, -90), tolerance = 1e-12)
  expect_equal(result$elevation[1:4], rep(0, 4), tolerance = 1e-12)
  expect_equal(result$elevation[5], 90, tolerance = 1e-12)

  center <- c(10, 20, 30)
  shifted <- cartesian_to_view_angles(
    x = center[1] + c(0, 1, 0, -1),
    y = center[2] + c(-1, 0, 1, 0),
    z = center[3] + c(0, 0, 0, 0),
    center = center
  )
  expect_equal(shifted$azimuth, c(0, 90, -180, -90), tolerance = 1e-12)
})


test_that("view_eye_face_on uses mean normals and places the long face-on range on y", {
  state <- get(".cv3d_face_on_state", envir = asNamespace("CV3D"))
  rm(list = ls(state, all.names = TRUE), envir = state)

  eye <- data.frame(
    x = c(-2, -1, 0, 1, 2, 0),
    y = c(0, 0.5, 1, 0.5, 0, -1),
    z = c(0, 0.1, 0, -0.1, 0, 0),
    norm.x = 0,
    norm.y = 0,
    norm.z = 1
  )

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 400, height = 400)
  on.exit({
    grDevices::dev.off()
    unlink(png_file)
  }, add = TRUE)

  view <- view_eye_face_on(eye, projection = "2D")

  expect_equal(view$view_direction, c(0, 0, 1), tolerance = 1e-12)
  expect_true(view$rotated_90)
  expect_gt(view$reference_height, view$reference_width)
  expect_equal(colMeans(view$projected_coordinates), c(x = 0, y = 0), tolerance = 1e-12)

  basis <- rbind(
    view$x_direction,
    view$y_direction,
    view$view_direction
  )
  expect_equal(unname(basis %*% t(basis)), diag(3), tolerance = 1e-12)
  expect_equal(det(basis), 1, tolerance = 1e-12)
})


test_that("view_eye_face_on centres comparison clouds and applies reference-sized translations", {
  state <- get(".cv3d_face_on_state", envir = asNamespace("CV3D"))
  rm(list = ls(state, all.names = TRUE), envir = state)

  eye <- data.frame(
    x = c(-2, -1, 0, 1, 2, 0),
    y = c(0, 0.5, 1, 0.5, 0, -1),
    z = c(0, 0.1, 0, -0.1, 0, 0),
    norm.x = 0,
    norm.y = 0,
    norm.z = 1
  )
  eye_shifted <- eye
  eye_shifted$x <- eye_shifted$x + 100
  eye_shifted$y <- eye_shifted$y - 50
  eye_shifted$z <- eye_shifted$z + 25

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 500, height = 400)
  on.exit({
    grDevices::dev.off()
    unlink(png_file)
  }, add = TRUE)

  first <- view_eye_face_on(eye, projection = "2D")

  second <- view_eye_face_on(
    eye_shifted,
    projection = "2D",
    add = TRUE,
    translate_x = 1.1,
    col = "red"
  )
  usr_second <- graphics::par("usr")

  expected_shift <- c(x = 1.1 * first$reference_width, y = 0)
  expect_equal(colMeans(second$projected_coordinates), expected_shift, tolerance = 1e-12)

  recentered_second <- sweep(
    second$projected_coordinates,
    2,
    expected_shift,
    "-"
  )
  expect_equal(recentered_second, first$projected_coordinates, tolerance = 1e-12)
  # The redrawn plotting region must contain both layers. With asp = 1, the
  # numeric x span need not itself grow because the device aspect ratio can
  # already provide sufficient horizontal space.
  all_x <- c(
    first$projected_coordinates[, "x"],
    second$projected_coordinates[, "x"]
  )
  all_y <- c(
    first$projected_coordinates[, "y"],
    second$projected_coordinates[, "y"]
  )
  expect_lte(usr_second[1], min(all_x))
  expect_gte(usr_second[2], max(all_x))
  expect_lte(usr_second[3], min(all_y))
  expect_gte(usr_second[4], max(all_y))
  expect_equal(
    second$translation_3d,
    1.1 * first$reference_width * first$x_direction,
    tolerance = 1e-12
  )
})


test_that("view_eye_face_on can reuse a returned view without surface normals", {
  state <- get(".cv3d_face_on_state", envir = asNamespace("CV3D"))
  rm(list = ls(state, all.names = TRUE), envir = state)

  eye <- data.frame(
    x = c(-2, -1, 0, 1, 2),
    y = c(0, 1, 1.5, 1, 0),
    z = c(0, 0, 0, 0, 0),
    norm.x = 0,
    norm.y = 0,
    norm.z = 1
  )

  png_file <- tempfile(fileext = ".png")
  grDevices::png(png_file, width = 400, height = 400)
  on.exit({
    grDevices::dev.off()
    unlink(png_file)
  }, add = TRUE)

  view <- view_eye_face_on(eye, projection = "2D")
  overlay <- eye[c(2, 4), c("x", "y", "z")]

  overlay_view <- view_eye_face_on(
    overlay,
    projection = "2D",
    view = view,
    add = FALSE,
    col = "red",
    pch = 16
  )

  expect_equal(overlay_view$view_direction, view$view_direction, tolerance = 1e-12)
  expect_equal(overlay_view$x_direction, view$x_direction, tolerance = 1e-12)
  expect_equal(overlay_view$y_direction, view$y_direction, tolerance = 1e-12)
  expect_equal(colMeans(overlay_view$projected_coordinates), c(x = 0, y = 0), tolerance = 1e-12)
})


test_that("view_eye_face_on uses the same translation convention in 3D", {
  skip_if_not_installed("rgl")

  state <- get(".cv3d_face_on_state", envir = asNamespace("CV3D"))
  rm(list = ls(state, all.names = TRUE), envir = state)

  eye <- data.frame(
    x = c(-2, -1, 0, 1, 2, 0),
    y = c(0, 0.5, 1, 0.5, 0, -1),
    z = c(0, 0.1, 0, -0.1, 0, 0),
    norm.x = 0,
    norm.y = 0,
    norm.z = 1
  )

  old_options <- options(rgl.useNULL = TRUE)
  on.exit(options(old_options), add = TRUE)
  on.exit(try(rgl::close3d(), silent = TRUE), add = TRUE)

  first <- view_eye_face_on(eye, projection = "3D")
  second <- view_eye_face_on(
    eye,
    projection = "3D",
    add = TRUE,
    translate_x = 1.1,
    translate_y = -0.5
  )

  expect_equal(
    second$translation_3d,
    1.1 * first$reference_width * first$x_direction -
      0.5 * first$reference_height * first$y_direction,
    tolerance = 1e-12
  )
})
