test_that("edge detection distinguishes a hexagonal lattice boundary", {
  rows <- -4:4
  pts <- do.call(rbind, lapply(rows, function(r) {
    cols <- -4:4
    data.frame(
      x = cols + 0.5 * (r %% 2),
      y = r * sqrt(3) / 2,
      z = 0.002 * (cols^2 + r^2)
    )
  }))
  pts$ID <- sprintf("F%05d", seq_len(nrow(pts)))
  pts <- pts[, c("ID", "x", "y", "z")]

  e90 <- detect_facet_edges(pts, gap_threshold_deg = 90)
  e100 <- detect_facet_edges(pts, gap_threshold_deg = 100)

  centre <- which.min(pts$x^2 + pts$y^2)
  expect_false(e90$is_edge_facet[centre])
  expect_true(sum(e90$is_edge_facet) >= sum(e100$is_edge_facet))
  expect_true(any(e90$is_edge_facet))
})

test_that("edge-aware neighbour links are reciprocal and bounded", {
  rows <- -4:4
  pts <- do.call(rbind, lapply(rows, function(r) {
    cols <- -4:4
    data.frame(
      x = cols + 0.5 * (r %% 2),
      y = r * sqrt(3) / 2,
      z = 0.002 * (cols^2 + r^2)
    )
  }))
  pts$ID <- sprintf("F%05d", seq_len(nrow(pts)))
  pts <- pts[, c("ID", "x", "y", "z")]

  g <- find_neighbours_edge_aware(pts, edge_gap_threshold_deg = 90)
  expect_true(all(g$number_of_neighbours <= 6))
  expect_true(all(c("is_edge_facet", "edge_angular_gap_deg", "neighbour_core_spacing_um") %in% names(g)))

  parse_ids <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return(character(0))
    trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  }
  for (i in seq_len(nrow(g))) {
    for (j_id in parse_ids(g$neighbours[i])) {
      j <- match(j_id, g$ID)
      expect_false(is.na(j))
      expect_true(g$ID[i] %in% parse_ids(g$neighbours[j]))
    }
  }

  centre <- which.min(pts$x^2 + pts$y^2)
  expect_equal(g$number_of_neighbours[centre], 6)
})
