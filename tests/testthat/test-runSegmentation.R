# setup
copykit_obj_50kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '50kb'
)
copykit_obj_50kb <- runVst(copykit_obj_50kb)

copykit_obj_100kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '100kb'
)
copykit_obj_100kb <- runVst(copykit_obj_100kb)

copykit_obj_175kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '175kb'
)
copykit_obj_175kb <- runVst(copykit_obj_175kb)

copykit_obj_200kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '200kb'
)
copykit_obj_200kb <- runVst(copykit_obj_200kb)

copykit_obj_250kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '250kb'
)
copykit_obj_250kb <- runVst(copykit_obj_250kb)

copykit_obj_500kb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '500kb'
)
copykit_obj_500kb <- runVst(copykit_obj_500kb)

copykit_obj_1Mb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '1Mb'
)
copykit_obj_1Mb <- runVst(copykit_obj_1Mb)

copykit_obj_2Mb <- mock_bincounts(
  ncells = 10,
  ncells_diploid = 5,
  position_gain = 1:50,
  position_del = 200:250,
  resolution = '2.5Mb'
)
copykit_obj_2Mb <- runVst(copykit_obj_2Mb)

# test
test_that("Testing CopyKit runSegmentation for different resolutions: ", {
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_50kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_100kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_175kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_200kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_250kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_500kb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_1Mb), "CopyKit")
  expect_s4_class(copykit_obj <- runSegmentation(copykit_obj_2Mb), "CopyKit")
})
