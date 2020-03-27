context("Gene selection")
genes = paste0("Gene_", seq_len(500))
id = sample(c(TRUE, FALSE), replace = TRUE, length(genes))
genesnull = genes[id]
genesnonnull = genes[!id]
genes2 = paste0("Gene_", seq_len(10))
id2 = sample(c(TRUE, FALSE), replace = TRUE, length(genes2))
genes2null = genes2[id2]
genes2nonnull = genes2[!id2]
batch = sample(LETTERS[1:3], 20, replace = TRUE)
group = sample(1:3, 20, replace = TRUE)
Design = configExperiment(batch.config = c(0.5, 0.3, 0.2), 
                          group.config = c(0.6, 0.4), tot.samples = 400, 
                          batch = batch, group = group)
test_that("Gene selection selects correct number of genes",{
    expect_identical(
      length(selectGenes(0.2, Design, 200, genesnull, genesnonnull)), 200
      )
    expect_identical(
      length(selectGenes(0, Design, 200, genesnull, genesnonnull)), 200
    )
    expect_identical(
      length(selectGenes(1, Design, 200, genesnull, genesnonnull)), 200
    )
    expect_identical(
      length(selectGenes(0.5, Design, 200, genes2null, genes2nonnull)), 200
    )
    expect_identical(
      length(selectGenes(0, Design, 200, genes2null, genes2nonnull)), 200
    )
    expect_identical(
      length(selectGenes(1, Design, 200, genes2null, genes2nonnull)), 200
    )
          })
