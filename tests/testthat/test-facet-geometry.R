parse_neighbour_ids <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}


test_that("find_neighbours returns valid bounded neighbour lists", {
  facets <- data.frame(
    ID = as.character(seq_len(8)),
    x = c(1, -1, 0, 0, 0, 0, 1, -1),
    y = c(0, 0, 1, -1, 0, 0, 1, -1),
    z = c(0, 0, 0, 0, 1, -1, 1, -1)
  )

  result <- find_neighbours(
    facets,
    k_local = 3,
    knn_search = 5,
    max_neighbours = 3
  )

  expect_equal(nrow(result), nrow(facets))
  expect_true(all(result$number_of_neighbours <= 3))

  for (i in seq_len(nrow(result))) {
    neighbours <- parse_neighbour_ids(result$neighbours[i])
    expect_false(result$ID[i] %in% neighbours)
    expect_true(all(neighbours %in% result$ID))
    expect_equal(length(neighbours), result$number_of_neighbours[i])

    for (neighbour_id in neighbours) {
      j <- match(neighbour_id, result$ID)
      reciprocal <- parse_neighbour_ids(result$neighbours[j])
      expect_true(result$ID[i] %in% reciprocal)
    }
  }
})


test_that("find_neighbours uses six as the default maximum and NULL disables the cap", {
  data("cv3d_example_facets", package = "CV3D")

  bounded <- find_neighbours(cv3d_example_facets)
  unbounded <- find_neighbours(cv3d_example_facets, max_neighbours = NULL)

  expect_true(all(bounded$number_of_neighbours <= 6))
  expect_true(all(unbounded$number_of_neighbours >= 0))

  for (i in seq_len(nrow(unbounded))) {
    neighbours <- parse_neighbour_ids(unbounded$neighbours[i])
    for (neighbour_id in neighbours) {
      j <- match(neighbour_id, unbounded$ID)
      reciprocal <- parse_neighbour_ids(unbounded$neighbours[j])
      expect_true(as.character(unbounded$ID[i]) %in% reciprocal)
    }
  }
})


test_that("calculate_facet_size uses Euclidean spacing and 50:50 smoothing", {
  facets <- data.frame(
    ID = c("A", "B", "C"),
    x = c(0, 1, 3),
    y = 0,
    z = 0,
    neighbours = c("B", "A; C", "B")
  )

  result <- calculate_facet_size(facets)

  expect_equal(result$facet_size, c(1, 1.5, 2), tolerance = 1e-12)
  expect_equal(result$n_neighbours_used, c(1L, 2L, 1L))
  expect_equal(result$facet_size_smoothed, c(1.25, 1.5, 1.75), tolerance = 1e-12)
})


test_that("facet normals are returned for a realistic facet patch", {
  data("cv3d_example_facets", package = "CV3D")
  neighbours <- find_neighbours(cv3d_example_facets)

  normals <- get_facet_normals(
    neighbours,
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )

  expect_true(all(c("ID", "norm.x", "norm.y", "norm.z") %in% names(normals)))
  expect_true(all(normals$ID %in% as.character(neighbours$ID)))
  finite_rows <- stats::complete.cases(normals[, c("norm.x", "norm.y", "norm.z")])
  expect_gt(sum(finite_rows), 0L)
})


testthat::test_that("get_facet_normals returns unit vectors", {
  data("cv3d_example_facets", package = "CV3D")
  neighbours <- find_neighbours(cv3d_example_facets)
  normals <- get_facet_normals(neighbours, cores = 1, verbose = FALSE)
  lengths <- sqrt(normals$norm.x^2 + normals$norm.y^2 + normals$norm.z^2)
  lengths <- lengths[is.finite(lengths)]
  testthat::expect_true(length(lengths) > 0)
  testthat::expect_equal(lengths, rep(1, length(lengths)), tolerance = 1e-10)
})


test_that("find_neighbours exposes a minimum-neighbour fallback", {
  expect_equal(formals(find_neighbours)$min_neighbours, 3)
  expect_error(
    find_neighbours(cv3d_example_facets, min_neighbours = 7, max_neighbours = 6),
    "cannot exceed"
  )
})


test_that("find_neighbours validates discrete and logical control parameters", {
  data("cv3d_example_facets", package = "CV3D")
  expect_error(find_neighbours(cv3d_example_facets, center = NA), "center")
  expect_error(find_neighbours(cv3d_example_facets, k_local = 2.5), "positive integer")
  expect_error(find_neighbours(cv3d_example_facets, knn_search = Inf), "positive integer")
  expect_error(find_neighbours(cv3d_example_facets, max_neighbours = 3.5), "positive integer")
})
