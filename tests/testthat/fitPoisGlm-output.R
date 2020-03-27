context("Fast fitting a Poisson glm")
y = rpois(5, 2)
x = rnorm(5)
offs = log(runif(5,1,2))
test_that("fitGlmPois does the same as glm",{
    expect_equivalent(fitPoisGlm(y, x, 2, offs)$coef,
                    glm(y~x+I(x^2), family = "poisson", offset = offs)$coef)
  expect_equivalent(fitPoisGlm(y, x, 3, offs)$coef,
                    glm(y~x+I(x^2)+I(x^3), family = "poisson", offset = offs)$coef)
          })
