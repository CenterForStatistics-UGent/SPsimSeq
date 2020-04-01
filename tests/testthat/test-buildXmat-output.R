context("Building the design matrix")
x = rnorm(5)
test_that("buildXmat does the same as model.matrix",{
    expect_equivalent(model.matrix(~I(x)+I(x^2)),
                    buildXmat(x, 3))
    expect_equivalent(model.matrix(~I(x)+I(x^2)+I(x^3)),
                    buildXmat(x, 4))
    expect_equivalent(model.matrix(~I(x)+I(x^2)+I(x^3)+I(x^4)),
                    buildXmat(x, 5))
          })
