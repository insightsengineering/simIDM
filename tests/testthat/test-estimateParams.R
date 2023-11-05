# prepareData ----

test_that("prepareData works as expected", {
  # Create simIDM data for the 4 possible different transition scenarios.
  colnames <- c(
    "id", "trt", "PFStime", "CensoredPFS", "PFSevent", "OStime",
    "CensoredOS", "OSevent", "recruitTime", "OStimeCal", "PFStimeCal"
  )
  patCens1 <- c(1, 1, 0.2, 1, 0, 0.2, 1, 0, 0.1, 0.3, 0.3)
  patCens2 <- c(2, 1, 0.7, 0, 1, 0.8, 1, 0, 1.2, 2, 1.9)
  pat13 <- c(3, 1, 0.4, 0, 1, 0.4, 0, 1, 0.4, 0.8, 0.8)
  pat123 <- c(4, 1, 0.1, 0, 1, 0.25, 0, 1, 0, 0.25, 0.1)
  df <- setNames(data.frame(rbind(patCens1, patCens2, pat13, pat123)), nm = colnames)
  actual <- prepareData(df)
  expect_snapshot(actual)
})
