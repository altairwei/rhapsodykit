test_that("test reassign_sample", {
  object <- list(
    sample = c(
      rep("sample-A-1", 5),
      rep("sample-A-2", 5),
      rep("sample-A-3", 5),
      rep("sample-B-1", 5),
      rep("sample-B-2", 5),
      rep("sample-B-3", 5),
      rep("sample-C-1", 5),
      rep("sample-C-2", 5),
      rep("sample-C-3", 5)
    )
  )

  new_obj <- reassign_sample(object, c(
    "sample-A-1" = "sample-A",
    "sample-A-2" = "sample-A",
    "sample-A-3" = "sample-A",
    "sample-B-1" = "sample-B",
    "sample-B-2" = "sample-B",
    "sample-B-3" = "sample-B",
    "sample-C-1" = "sample-C",
    "sample-C-2" = "sample-C",
    "sample-C-3" = "sample-C"
  ))

  expect_equal(new_obj$sample, c(
      rep("sample-A", 15),
      rep("sample-B", 15),
      rep("sample-C", 15)
    )
  )

  new_obj2 <- reassign_sample(object, c(
    "sample-A-1" = "sample-A",
    "sample-A-2" = "sample-A",
    "sample-A-3" = "sample-A"
  ))

  expect_equal(new_obj2$sample, c(
      rep("sample-A", 15),
      rep("sample-B-1", 5),
      rep("sample-B-2", 5),
      rep("sample-B-3", 5),
      rep("sample-C-1", 5),
      rep("sample-C-2", 5),
      rep("sample-C-3", 5)
    )
  )

  new_obj3 <- reassign_sample(object, c(
    "sample-A-1" = "sample-A",
    "sample-B-2" = "sample-B",
    "sample-C-3" = "sample-C"
  ))

  expect_equal(new_obj3$sample, c(
      rep("sample-A", 5),
      rep("sample-A-2", 5),
      rep("sample-A-3", 5),
      rep("sample-B-1", 5),
      rep("sample-B", 5),
      rep("sample-B-3", 5),
      rep("sample-C-1", 5),
      rep("sample-C-2", 5),
      rep("sample-C", 5)
    )
  )
})
