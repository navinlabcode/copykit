set.seed(1000)
copykit_obj <- copykit_example_filtered()[,sample(200)]
copykit_obj <- findClusters(copykit_obj)
copykit_obj <- runPca(copykit_obj)

test_that("Testing reduced embedding plots: ", {
  # testing UMAP
    expect_s3_class(p <- plotUmap(copykit_obj), "ggplot")
    expect_s3_class(p <- plotUmap(copykit_obj, label = "subclones"), "ggplot")

    # testing PCA
    expect_s3_class(p <- plotPca(copykit_obj), 'ggplot')
    expect_s3_class(p <- plotPca(copykit_obj, label = 'subclones'), 'ggplot')
})
