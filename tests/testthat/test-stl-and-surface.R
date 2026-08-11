test_that("STL_triangles parses a simple ASCII triangle", {
  stl_file <- tempfile(fileext = ".stl")
  on.exit(unlink(stl_file), add = TRUE)

  writeLines(c(
    "solid example",
    "facet normal 0 0 1",
    "outer loop",
    "vertex 0 0 0",
    "vertex 1 0 0",
    "vertex 0 1 0",
    "endloop",
    "endfacet",
    "endsolid example"
  ), stl_file)

  triangles <- STL_triangles(stl_file, plot_results = FALSE, verbose = FALSE)

  expect_s3_class(triangles, "tbl_df")
  expect_equal(nrow(triangles), 1L)
  expect_equal(names(triangles), c("ID", "x", "y", "z", "norm.x", "norm.y", "norm.z"))
  expect_equal(triangles$ID, 1)
  expect_equal(c(triangles$x, triangles$y, triangles$z), c(1 / 3, 1 / 3, 0), tolerance = 1e-12)
  expect_equal(c(triangles$norm.x, triangles$norm.y, triangles$norm.z), c(0, 0, 1))
})


test_that("calculate_local_heights returns zero heights for a planar surface", {
  surface <- expand.grid(x = -1:1, y = -1:1)
  surface$z <- 0
  surface$norm.x <- 0
  surface$norm.y <- 0
  surface$norm.z <- 1

  result <- calculate_local_heights(
    surface,
    neighbourhood_radius = 1,
    cores = 1,
    verbose = FALSE
  )

  expect_equal(nrow(result), nrow(surface))
  expect_true(all(c(
    "local_height", "local_height_col", "local_height_filtered_col",
    "local_height_contrast", "local_height_contrast_col"
  ) %in% names(result)))
  expect_equal(result$local_height, rep(0, nrow(surface)), tolerance = 1e-12)
  expect_true(all(result$local_height_contrast >= 0 & result$local_height_contrast <= 1))
})


test_that("normalize_local_heights produces bounded local values", {
  surface <- expand.grid(x = -1:1, y = -1:1)
  surface$z <- 0
  surface$local_height <- seq_len(nrow(surface))

  result <- normalize_local_heights(
    surface,
    neighbourhood_radius = 1.5,
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )

  expect_equal(nrow(result), nrow(surface))
  expect_true(all(c(
    "local_height_norm", "local_height_norm_col",
    "local_height_norm_contrast", "local_height_norm_contrast_col"
  ) %in% names(result)))
  finite_norm <- is.finite(result$local_height_norm)
  expect_true(any(finite_norm))
  expect_true(all(result$local_height_norm[finite_norm] >= 0 & result$local_height_norm[finite_norm] <= 1))
  expect_true(all(result$local_height_norm_contrast[finite_norm] >= 0 & result$local_height_norm_contrast[finite_norm] <= 1))
  expect_equal(result$local_height_norm_contrast[finite_norm],
               (10^result$local_height_norm[finite_norm] - 1) / 9, tolerance = 1e-12)
})


test_that("spherical local-height neighbourhoods are rotation invariant", {
  surface <- expand.grid(x = -1:1, y = -1:1)
  surface$z <- 0.15 * surface$x^2 + 0.08 * surface$y
  surface$norm.x <- 0
  surface$norm.y <- 0
  surface$norm.z <- 1

  angle <- pi / 4
  rotation <- matrix(c(
    cos(angle), -sin(angle), 0,
    sin(angle),  cos(angle), 0,
    0,           0,          1
  ), nrow = 3, byrow = TRUE)

  rotated <- surface
  rotated[, c("x", "y", "z")] <- as.matrix(surface[, c("x", "y", "z")]) %*% t(rotation)
  rotated[, c("norm.x", "norm.y", "norm.z")] <-
    as.matrix(surface[, c("norm.x", "norm.y", "norm.z")]) %*% t(rotation)

  original_heights <- calculate_local_heights(
    surface, neighbourhood_radius = 1.1, cores = 1, verbose = FALSE
  )
  rotated_heights <- calculate_local_heights(
    rotated, neighbourhood_radius = 1.1, cores = 1, verbose = FALSE
  )

  expect_equal(
    rotated_heights$local_height,
    original_heights$local_height,
    tolerance = 1e-10
  )

  original_norm <- normalize_local_heights(
    original_heights,
    neighbourhood_radius = 1.5,
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )
  rotated_norm <- normalize_local_heights(
    rotated_heights,
    neighbourhood_radius = 1.5,
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )

  expect_equal(
    rotated_norm$local_height_norm,
    original_norm$local_height_norm,
    tolerance = 1e-10
  )
})


test_that("normalization leaves unsupported points missing and exposes quantiles", {
  surface <- data.frame(
    x = c(0, 0.1, 100),
    y = c(0, 0, 0),
    z = c(0, 0, 0),
    local_height = c(0, 1, 5)
  )
  result <- normalize_local_heights(
    surface,
    neighbourhood_radius = 0.2,
    lower_quantile = 0.10,
    upper_quantile = 0.90,
    cores = 1,
    verbose = FALSE
  )
  expect_true(is.na(result$local_height_norm[3]))
  expect_true(is.na(result$local_height_norm_contrast[3]))
  expect_equal(formals(normalize_local_heights)$lower_quantile, 0.10)
  expect_equal(formals(normalize_local_heights)$upper_quantile, 0.90)
})
