test_that("Function creates a graph that shows the number of global and regional observations in each grid", {
  load("data/JoinnedDF.rda")
  expect_is(narccapCOUNT(JoinnedDF, 1986, 8), "ggplot")
})
