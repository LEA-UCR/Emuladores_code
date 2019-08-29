context("test-ShowData")

test_that("ShowDataFunctionWorks", {
  expect_is(ShowData(), "data.frame")
})
