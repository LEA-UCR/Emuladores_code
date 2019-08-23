test_that("Function creates a map", {
  load("data/sic Regional.rda")
  expect_is(narccapMAP(`sic Regional`, "sic", 1985, 1), "ggplot")
})
