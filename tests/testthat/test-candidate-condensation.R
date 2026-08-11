test_that("facet condensation identifies separated candidate groups", {
  points <- data.frame(
    x = c(-0.2, 0, 0.2, 0, 0, 9.8, 10, 10.2, 10, 10),
    y = c(0, 0, 0, -0.2, 0.2, 0, 0, 0, -0.2, 0.2),
    z = 0,
    height_value = c(1, 4, 1, 2, 2, 1, 5, 1, 2, 2)
  )

  result <- find_facet_candidates_condensed(
    points,
    neighbour_radius = 0.5,
    merge_radius = 0.3,
    min_cluster_size = 3,
    cores = 1,
    verbose = FALSE
  )

  expect_equal(nrow(result), 2L)
  expect_true(all(result$selected_source_index %in% seq_len(nrow(points))))

  for (i in seq_len(nrow(result))) {
    source_row <- points[result$selected_source_index[i], c("x", "y", "z")]
    expect_equal(
      as.numeric(result[i, c("x", "y", "z")]),
      as.numeric(source_row),
      tolerance = 1e-12
    )
  }
})


test_that("facet condensation has consistent detailed empty output", {
  empty <- data.frame(x = numeric(), y = numeric(), z = numeric(), height_value = numeric())

  result <- find_facet_candidates_condensed(
    empty,
    neighbour_radius = 1,
    return_details = TRUE,
    cores = 1,
    verbose = FALSE
  )

  expect_named(result, c("candidates", "membership", "parameters"))
  expect_equal(nrow(result$candidates), 0L)
  expect_equal(nrow(result$membership), 0L)
  expect_true(isTRUE(result$parameters$converged))
})


test_that("facet condensation does not merge transitive chains beyond merge_radius", {
  points <- data.frame(
    x = c(0, 0.2, 0.4),
    y = 0,
    z = 0,
    height_value = 1
  )

  result <- find_facet_candidates_condensed(
    points,
    neighbour_radius = 1e-6,
    merge_radius = 0.3,
    max_iterations = 1,
    min_cluster_size = 1,
    cores = 1,
    verbose = FALSE,
    return_details = TRUE
  )

  memberships <- split(result$membership$x, result$membership$cluster_id)
  max_spans <- vapply(memberships, function(x) max(x) - min(x), numeric(1))
  expect_true(all(max_spans <= 0.3 + 1e-12))
  expect_gt(length(memberships), 1)
})


test_that("non-negative candidate heights retain their absolute weighting scale", {
  points <- data.frame(
    x = c(0, 1, 2),
    y = 0,
    z = 0,
    height_value = c(0.6, 0.7, 0.9)
  )

  result <- find_facet_candidates_condensed(
    points,
    neighbour_radius = 3,
    merge_radius = 3,
    weight_exponent = 2,
    max_iterations = 1,
    step_size = 1,
    min_cluster_size = 1,
    cores = 1,
    verbose = FALSE
  )

  expected_mode_x <- sum(points$x * points$height_value^2) /
    sum(points$height_value^2)

  expect_equal(result$mode_x, expected_mode_x, tolerance = 1e-12)
})
