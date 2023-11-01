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

  expect_data_frame(actual, nrows = 10, ncols = 7)
  expect_equal(names(actual), c("id", "from", "to", "trans", "entry", "exit", "status"))
  expect_equal(as.numeric(actual[1, ]), c(1, 1, 2, 1, 0, 0.2, 0))
  expect_equal(as.numeric(actual[2, ]), c(1, 1, 3, 2, 0, 0.2, 0))
  expect_equal(as.numeric(actual[3, ]), c(2, 1, 2, 1, 0, 0.7, 1))
  expect_equal(as.numeric(actual[4, ]), c(2, 1, 3, 2, 0, 0.7, 0))
  expect_equal(as.numeric(actual[5, ]), c(2, 2, 3, 3, 0.7, 0.8, 0))
  expect_equal(as.numeric(actual[6, ]), c(3, 1, 2, 1, 0, 0.4, 0))
  expect_equal(as.numeric(actual[7, ]), c(3, 1, 3, 2, 0, 0.4, 1))
  expect_equal(as.numeric(actual[8, ]), c(4, 1, 2, 1, 0, 0.1, 1))
  expect_equal(as.numeric(actual[9, ]), c(4, 1, 3, 2, 0, 0.1, 0))
  expect_equal(as.numeric(actual[10, ]), c(4, 2, 3, 3, 0.1, 0.25, 1))
})
