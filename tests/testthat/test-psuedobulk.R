test_that("test make_pseudo_bulk", {
  mat_list <- list(
    sample1 = Matrix::Matrix(
      c(
        1, 2, 3, # 2
        4, 5, 6, # 5
        7, 8, 9, # 8
        2, 4, 6  # 4
      ),
      dimnames = list(
        c("gene1", "gene2", "gene3", "gene4"), # rownames
        c("cell1", "cell2", "cell3") # colnames
      ),
      byrow = TRUE,
      nrow = 4, ncol = 3,
    ),
    sample2 = Matrix::Matrix(
      c(
        6, 2, 1, # 3
        4, 5, 6, # 5
        7, 8, 3, # 6
        2, 4, 3, # 3
        3, 7, 2, # 4
        3, 4, 2  # 3
      ),
      dimnames = list(
        c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6"), # rownames
        c("cell1", "cell2", "cell3") # colnames
      ),
      byrow = TRUE,
      nrow = 6, ncol = 3,
    ),
    sample3 = Matrix::Matrix(
      c(
        6, 2, 1, # 3
        4, 8, 6  # 6
      ),
      dimnames = list(
        c("gene1", "gene9"), # rownames
        c("cell1", "cell2", "cell3") # colnames
      ),
      byrow = TRUE,
      nrow = 2, ncol = 3,
    )
  )

  df <- .make_pseudo_bulk(mat_list, "avg")
  expect_equal(
    names(df),
    c("sample1", "sample2", "sample3"))
  expect_equal(
    rownames(df),
    c("gene1", "gene2", "gene3", "gene4", "gene5", "gene6", "gene9"))
  expect_equal(
    df[c("gene5", "gene6", "gene9"), "sample1"],
    c(0, 0, 0)
  )
  expect_equal(
    df["gene9", "sample2"],
    0
  )
  expect_equal(
    df[c("gene2", "gene3", "gene4", "gene5", "gene6"), "sample3"],
    c(0, 0, 0, 0, 0)
  )
})
