copykit_obj <- copykit_example_filtered()
copykit_obj <- findClusters(copykit_obj)

test_that("Testing reduced embedding plots: ", {
  expect_s3_class(p <- plotUmap(copykit_obj), "ggplot")
  expect_s3_class(p <- plotUmap(copykit_obj, label = 'subclones'), "ggplot")
})
