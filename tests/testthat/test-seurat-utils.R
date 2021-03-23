test_that("test assign_stim", {
  object <- list(
    stim = c(
      rep("stim-A-1", 5),
      rep("stim-A-2", 5),
      rep("stim-A-3", 5),
      rep("stim-B-1", 5),
      rep("stim-B-2", 5),
      rep("stim-B-3", 5),
      rep("stim-C-1", 5),
      rep("stim-C-2", 5),
      rep("stim-C-3", 5)
    )
  )

  new_obj <- assign_stim(object, c(
    "stim-A-1" = "stim-A",
    "stim-A-2" = "stim-A",
    "stim-A-3" = "stim-A",
    "stim-B-1" = "stim-B",
    "stim-B-2" = "stim-B",
    "stim-B-3" = "stim-B",
    "stim-C-1" = "stim-C",
    "stim-C-2" = "stim-C",
    "stim-C-3" = "stim-C"
  ))

  expect_equal(new_obj$stim, c(
      rep("stim-A", 15),
      rep("stim-B", 15),
      rep("stim-C", 15)
    )
  )

  new_obj2 <- assign_stim(object, c(
    "stim-A-1" = "stim-A",
    "stim-A-2" = "stim-A",
    "stim-A-3" = "stim-A"
  ))

  expect_equal(new_obj2$stim, c(
      rep("stim-A", 15),
      rep("stim-B-1", 5),
      rep("stim-B-2", 5),
      rep("stim-B-3", 5),
      rep("stim-C-1", 5),
      rep("stim-C-2", 5),
      rep("stim-C-3", 5)
    )
  )

  new_obj3 <- assign_stim(object, c(
    "stim-A-1" = "stim-A",
    "stim-B-2" = "stim-B",
    "stim-C-3" = "stim-C"
  ))

  expect_equal(new_obj3$stim, c(
      rep("stim-A", 5),
      rep("stim-A-2", 5),
      rep("stim-A-3", 5),
      rep("stim-B-1", 5),
      rep("stim-B", 5),
      rep("stim-B-3", 5),
      rep("stim-C-1", 5),
      rep("stim-C-2", 5),
      rep("stim-C", 5)
    )
  )
})
