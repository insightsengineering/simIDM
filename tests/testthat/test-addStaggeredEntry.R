# addStaggeredEntry ----
test_that("addStaggeredEntry works as expected if no staggered study entry is present", {
  SimData <- data.frame(
    id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2), entry = c(0, 3, 0, 0),
    exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
  )
  actual <- addStaggeredEntry(SimData, N = 3, accrualParam = "intensity", accrualValue = 0)
  expected <- cbind(SimData, entryAct = SimData$entry, exitAct = SimData$exit, censAct = SimData$censTime)
  expected[, c("censTime")] <- list(NULL)
  expect_equal(actual, expected)


  actual2 <- addStaggeredEntry(SimData, N = 3, accrualParam = "time", accrualValue = 0)
  expected2 <- cbind(SimData, entryAct = SimData$entry, exitAct = SimData$exit, censAct = SimData$censTime)
  expected2[, c("censTime")] <- list(NULL)
  expect_equal(actual2, expected2)
})


test_that("addStaggeredEntry works as expected if staggered study entry is specified by accrual parameter time", {
  SimData <- data.frame(
    id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2), entry = c(0, 3, 0, 0),
    exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
  )
  actual <- addStaggeredEntry(SimData, N = 3, accrualParam = "time", accrualValue = 6)

  # accrual times have to be between 0 and 6.
  expect_true(all((actual$entryAct - actual$entry) <= 6))
  expect_true(all((actual$entryAct - actual$entry) >= 0))
  expect_true(all(actual$entryAct > actual$entry))

  # the actual exit time is the individual exit time + the actual entry time at study scale.
  expect_equal(actual$exitAct - (actual$entryAct - actual$entry), actual$exit)
})



test_that("addStaggeredEntry works as expected if staggered study entry is specified by accrual parameter intensity", {
  SimData <- data.frame(
    id = c(1, 1, 2, 3), from = c(0, 1, 0, 0), to = c(1, 2, "cens", 2), entry = c(0, 3, 0, 0),
    exit = c(3, 5.3, 5.6, 7.2), censTime = c(6.8, 6.8, 5.6, 9.4)
  )
  actual <- addStaggeredEntry(SimData, N = 3, accrualParam = "intensity", accrualValue = 8)

  # Accrual times have to be between 0 and N/accrualValue.
  expValue <- 3 / 8
  expect_true(all((actual$entryAct - actual$entry) <= expValue))
  expect_true(all((actual$entryAct - actual$entry) >= 0))
  expect_true(all(actual$entryAct > actual$entry))

  # The actual exit time is the individual exit time + the actual entry time at study scale.
  expect_equal(actual$exitAct - (actual$entryAct - actual$entry), actual$exit)
})
