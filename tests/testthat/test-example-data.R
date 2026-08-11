test_that("bundled example datasets have the documented compact structure", {
  data("cv3d_example_surface", package = "CV3D")
  data("cv3d_example_thresholded", package = "CV3D")
  data("cv3d_example_facets", package = "CV3D")

  expect_equal(nrow(cv3d_example_surface), 600L)
  expect_equal(names(cv3d_example_surface), c("ID", "x", "y", "z", "norm.x", "norm.y", "norm.z"))
  expect_equal(cv3d_example_surface$ID, seq_len(nrow(cv3d_example_surface)))

  expect_equal(nrow(cv3d_example_thresholded), 74L)
  expect_equal(names(cv3d_example_thresholded), c("source_index", "x", "y", "z", "height_value"))
  expect_equal(cv3d_example_thresholded$source_index, seq_len(nrow(cv3d_example_thresholded)))

  expect_equal(nrow(cv3d_example_facets), 61L)
  expect_equal(names(cv3d_example_facets), c("ID", "x", "y", "z"))
  expect_equal(cv3d_example_facets$ID, seq_len(nrow(cv3d_example_facets)))
})
