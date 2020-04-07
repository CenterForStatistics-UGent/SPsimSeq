context("SPsimSeq input")
data("zhang.data.sub")
zhang.counts <- zhang.data.sub$counts
MYCN.status  <- zhang.data.sub$MYCN.status
test_that("SPsimSeq throws errors for wrong input types", {
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 20,
                        pDE = 0.1, lfc.thrld = 0.5, result.format = "data.frame"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 1,
                        group.config = c(0.25, 0.5), tot.samples = 20,
                        pDE = 0.1, lfc.thrld = 0.5, result.format = "list"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 0.9,
                        group.config = c(0.5, 0.5), tot.samples = 20,
                        pDE = 0.1, lfc.thrld = 0.5, result.format = "list"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 20,
                        pDE = 20, lfc.thrld = 0.5, result.format = "list"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 1,
                        group.config = c(0.1, 0.1, 0.5, 0.3), tot.samples = 20,
                        pDE = 0.2, lfc.thrld = 0.5, result.format = "list"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, group.config = 1,
                        batch.config = c(0.1, 0.1, 0.5, 0.3), tot.samples = 20,
                        pDE = 0.2, lfc.thrld = 0.5, result.format = "list"))
  expect_error(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 2000, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 20,
                        pDE = 0.2, lfc.thrld = 0.5, result.format = "list",
                        lib.size.params = c(10, 5)))
})

test_that("SPsimSeq throws warnings when log-fold change threshold too high", {
  expect_warning(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                        group = MYCN.status, n.genes = 20, batch.config = 1,
                        group.config = c(0.5, 0.5), tot.samples = 8,
                        pDE = 0.1, lfc.thrld = 50, result.format = "list"))
})