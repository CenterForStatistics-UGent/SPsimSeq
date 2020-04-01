context("Configuring the experiment")
batch = sample(LETTERS[1:3], 20, replace = TRUE)
group = sample(1:3, 20, replace = TRUE)
test_that("configExperiment yields desired output",{
    expect_equivalent(configExperiment(batch.config=1, group.config=1, tot.samples=100, 
                                            batch = batch, group = group)$exprmt.config,
    100)
    expect_equivalent(configExperiment(batch.config=c(0.5,0.3,0.2), group.config=c(0.6,0.4), tot.samples = 400, 
                                             batch = batch, group = group)$exprmt.config,
                            matrix(c(120, 72, 48, 80, 48, 32), nrow = 3, ncol = 2))
          })
