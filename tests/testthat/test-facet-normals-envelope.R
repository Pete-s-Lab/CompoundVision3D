test_that("envelope normals are finite unit vectors and outward on an icosahedron", {
  phi <- (1 + sqrt(5)) / 2
  xyz <- rbind(
    c(0,  1,  phi), c(0, -1,  phi), c(0,  1, -phi), c(0, -1, -phi),
    c( 1,  phi, 0), c(-1,  phi, 0), c( 1, -phi, 0), c(-1, -phi, 0),
    c( phi, 0,  1), c( phi, 0, -1), c(-phi, 0,  1), c(-phi, 0, -1)
  )
  xyz <- xyz / sqrt(rowSums(xyz^2))
  ids <- as.character(seq_len(nrow(xyz)))

  d <- as.matrix(stats::dist(xyz))
  diag(d) <- Inf
  neighbour_idx <- lapply(seq_len(nrow(xyz)), function(i) order(d[i, ])[1:5])

  facets <- tibble::tibble(
    ID = ids,
    x = xyz[, 1],
    y = xyz[, 2],
    z = xyz[, 3],
    neighbours = vapply(
      neighbour_idx,
      function(ii) paste(ids[ii], collapse = "; "),
      character(1)
    )
  )

  out <- get_facet_normals_envelope(facets, envelope_factor = 1.25)
  expect_equal(nrow(out), nrow(facets))
  expect_true(all(stats::complete.cases(out[, c("norm.x", "norm.y", "norm.z")])))

  nrm <- as.matrix(out[, c("norm.x", "norm.y", "norm.z")])
  expect_equal(sqrt(rowSums(nrm^2)), rep(1, nrow(nrm)), tolerance = 1e-10)
  expect_true(all(rowSums(nrm * xyz) > 0))
  expect_true(all(out$normal_support_face_count > 0))
  expect_true(all(out$normal_method == "envelope"))
  expect_equal(unique(out$normal_envelope_factor), 1.25)
})

test_that("envelope normals are translation invariant", {
  phi <- (1 + sqrt(5)) / 2
  xyz <- rbind(
    c(0,  1,  phi), c(0, -1,  phi), c(0,  1, -phi), c(0, -1, -phi),
    c( 1,  phi, 0), c(-1,  phi, 0), c( 1, -phi, 0), c(-1, -phi, 0),
    c( phi, 0,  1), c( phi, 0, -1), c(-phi, 0,  1), c(-phi, 0, -1)
  )
  xyz <- xyz / sqrt(rowSums(xyz^2))
  ids <- as.character(seq_len(nrow(xyz)))
  d <- as.matrix(stats::dist(xyz))
  diag(d) <- Inf
  neighbour_idx <- lapply(seq_len(nrow(xyz)), function(i) order(d[i, ])[1:5])
  neighbours <- vapply(neighbour_idx, function(ii) paste(ids[ii], collapse = "; "), character(1))

  a <- tibble::tibble(ID = ids, x = xyz[, 1], y = xyz[, 2], z = xyz[, 3], neighbours = neighbours)
  b <- a
  b$x <- b$x + 123.4
  b$y <- b$y - 51.2
  b$z <- b$z + 9.8

  na <- get_facet_normals_envelope(a, envelope_factor = 1.25)
  nb <- get_facet_normals_envelope(b, envelope_factor = 1.25)

  expect_equal(
    as.matrix(na[, c("norm.x", "norm.y", "norm.z")]),
    as.matrix(nb[, c("norm.x", "norm.y", "norm.z")]),
    tolerance = 1e-9
  )
})

test_that("supported envelope factors return valid normals", {
  phi <- (1 + sqrt(5)) / 2
  xyz <- rbind(
    c(0,  1,  phi), c(0, -1,  phi), c(0,  1, -phi), c(0, -1, -phi),
    c( 1,  phi, 0), c(-1,  phi, 0), c( 1, -phi, 0), c(-1, -phi, 0),
    c( phi, 0,  1), c( phi, 0, -1), c(-phi, 0,  1), c(-phi, 0, -1)
  )
  xyz <- xyz / sqrt(rowSums(xyz^2))
  ids <- as.character(seq_len(nrow(xyz)))
  d <- as.matrix(stats::dist(xyz)); diag(d) <- Inf
  neighbour_idx <- lapply(seq_len(nrow(xyz)), function(i) order(d[i, ])[1:5])
  facets <- tibble::tibble(
    ID = ids, x = xyz[, 1], y = xyz[, 2], z = xyz[, 3],
    neighbours = vapply(neighbour_idx, function(ii) paste(ids[ii], collapse = "; "), character(1))
  )

  for (factor in c(1, 1.25, 1.5, 2)) {
    out <- get_facet_normals_envelope(facets, envelope_factor = factor)
    expect_true(all(is.finite(out$norm.x)))
    expect_equal(unique(out$normal_envelope_factor), factor)
  }
  expect_error(get_facet_normals_envelope(facets, envelope_factor = 0), "greater than zero")
})

test_that("envelope normals use facet positions directly and default to 1.25x", {
  phi <- (1 + sqrt(5)) / 2
  xyz <- rbind(
    c(0,  1,  phi), c(0, -1,  phi), c(0,  1, -phi), c(0, -1, -phi),
    c( 1,  phi, 0), c(-1,  phi, 0), c( 1, -phi, 0), c(-1, -phi, 0),
    c( phi, 0,  1), c( phi, 0, -1), c(-phi, 0,  1), c(-phi, 0, -1)
  )
  xyz <- xyz / sqrt(rowSums(xyz^2))
  facets <- tibble::tibble(
    ID = as.character(seq_len(nrow(xyz))),
    x = xyz[, 1], y = xyz[, 2], z = xyz[, 3]
  )

  out <- get_facet_normals_envelope(facets)
  expect_equal(nrow(out), nrow(facets))
  expect_equal(unique(out$normal_envelope_factor), 1.25)
  expect_true(all(out$normal_method == "envelope"))
})
