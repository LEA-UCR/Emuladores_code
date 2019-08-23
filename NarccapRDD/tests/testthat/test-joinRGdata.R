test_that("Function joins data frames", {
  load("data/PRECL Global.rda")
  load("data/sic Regional.rda")
  expect_is(joinRGdata(`PRECL Global`, `sic Regional`), "data.frame")
})

