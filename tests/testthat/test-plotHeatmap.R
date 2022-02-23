# setup
copykit_obj <- copykit_example_filtered()
set.seed(1000)
copykit_obj <- copykit_obj[, sample(70)]
copykit_obj <- findClusters(copykit_obj)
copykit_obj <- calcConsensus(copykit_obj)
copykit_obj <- runConsensusPhylo(copykit_obj)

#tests

test_that("Testing plotting heatmap: ", {
  expect_s4_class(ht <- plotHeatmap(copykit_obj), "Heatmap")

  # checking visual parameters
  expect_s4_class(ht <-
                    plotHeatmap(copykit_obj, label = "subclones"),
                  "Heatmap")

  expect_s4_class(ht <- plotHeatmap(copykit_obj, genes = c("MYC",
                                                           "TP53")),
                  "Heatmap")
})
