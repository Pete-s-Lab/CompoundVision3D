test_that("vector_angle returns geometrically correct angles", {
  expect_equal(CV3D:::vector_angle(c(1, 0, 0), c(0, 1, 0)), pi / 2, tolerance = 1e-12)
  expect_equal(CV3D:::vector_angle(c(1, 0, 0), c(0, 1, 0), unit = "degrees"), 90, tolerance = 1e-12)
  expect_true(is.na(CV3D:::vector_angle(c(0, 0, 0), c(1, 0, 0))))
})


test_that("calculate_normal returns unit and degenerate triangle normals", {
  normal <- CV3D:::calculate_normal(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  expect_equal(normal, c(0, 0, 1), tolerance = 1e-12)

  degenerate <- CV3D:::calculate_normal(c(0, 0, 0), c(1, 0, 0), c(2, 0, 0))
  expect_equal(degenerate, c(0, 0, 0), tolerance = 1e-12)
})


test_that("get_optic_properties applies lattice-specific Snyder metrics and local anatomical acuity", {
  theta <- 10 * pi / 180
  facets <- data.frame(
    ID = c("A", "B"),
    neighbours = c("B", "A"),
    norm.x = c(1, cos(theta)),
    norm.y = c(0, sin(theta)),
    norm.z = c(0, 0),
    facet_size_smoothed = c(2, 2)
  )

  hex <- get_optic_properties(
    facets,
    lattice = "hexagonal",
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )
  square <- get_optic_properties(
    facets,
    lattice = "square",
    cores = 1,
    plot_results = FALSE,
    verbose = FALSE
  )

  expect_equal(hex$interfacet_angle_deg, c(10, 10), tolerance = 1e-10)
  expect_equal(hex$interfacet_angle_rad, rep(theta, 2), tolerance = 1e-12)
  expect_equal(hex$sampling_lattice, c("hexagonal", "hexagonal"))
  expect_equal(square$sampling_lattice, c("square", "square"))

  expect_equal(
    hex$eye_parameter,
    rep(2 * theta * sqrt(3) / 2, 2),
    tolerance = 1e-12
  )
  expect_equal(
    square$eye_parameter,
    rep(2 * theta, 2),
    tolerance = 1e-12
  )
  expect_equal(
    hex$sampling_frequency_rad,
    rep(1 / (sqrt(3) * theta), 2),
    tolerance = 1e-12
  )
  expect_equal(
    square$sampling_frequency_rad,
    rep(1 / (2 * theta), 2),
    tolerance = 1e-12
  )

  # Local anatomical acuity uses 1 / (2 * delta_phi_deg), following the interommatidial-angle relationship used by Feller et al. (2020).
  # and is not lattice-corrected.
  expect_equal(hex$acuity_cpd, c(0.05, 0.05), tolerance = 1e-12)
  expect_equal(square$acuity_cpd, c(0.05, 0.05), tolerance = 1e-12)
})


test_that("get_optic_properties validates the lattice argument", {
  facets <- data.frame(
    ID = c("A", "B"),
    neighbours = c("B", "A"),
    norm.x = c(1, 1),
    norm.y = c(0, 0),
    norm.z = c(0, 0),
    facet_size_smoothed = c(1, 1)
  )

  expect_error(
    get_optic_properties(facets, lattice = "triangular"),
    "should be one of"
  )
})
