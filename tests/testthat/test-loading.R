test_that("package loads and exported functions are visible", {
  # Estimators
  expect_true(is.function(sdr))
  expect_true(is.function(itmle))
  expect_true(is.function(qreg))

  # Building blocks
  expect_true(is.function(density_ratio))
  expect_true(is.function(absorb_rule))
  expect_true(is.function(CausalState:::expand_to_horizon))

  # At least one SL wrapper
  expect_true(is.function(SL.tgt.glm))
  expect_true(is.function(SL.tmle_glm))
})
