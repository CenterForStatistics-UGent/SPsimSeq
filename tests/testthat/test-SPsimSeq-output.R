context("Output of the SPsimSeq function")
data("zhang.data.sub")
zhang.counts <- zhang.data.sub$counts
MYCN.status  <- zhang.data.sub$MYCN.status
test_that("SPsimSeq returns correct output types", {
  expect_type(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                          group = MYCN.status, n.genes = 10, batch.config = 1,
                          group.config = c(0.5, 0.5), tot.samples = 4,
                          pDE = 0.1, lfc.thrld = 0.5, result.format = "list"),
              "list")
  expect_s4_class(SPsimSeq(n.sim = 1, s.data = zhang.counts,
                       group = MYCN.status, n.genes = 10, batch.config = 1,
                       group.config = c(0.5, 0.5), tot.samples = 4,
                       pDE = 0.1, lfc.thrld = 0.5, result.format = "SCE")[[1]],
              "SingleCellExperiment")
})
