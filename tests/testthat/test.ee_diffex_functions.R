################################################################################

context("Unit tests for ee_diffex_functions.R")

################################################################################

#' @importFrom   tibble        data_frame
#'
test_that(".get_consistent_diffexes", {
  expect_error(
    object = .get_consistent_diffexes(
      dframe       = tibble::data_frame(gene   = letters[1:3],
                                        diffex = c(-1, 1, 0)
                                        ),
      grouping_col = "NOT PRESENT",
      diffex_col   = "diffex"
    ),
    info = "Grouping column is absent"
  )

  expect_error(
    object = .get_consistent_diffexes(
      dframe       = tibble::data_frame(gene.id = letters[1:3],
                                        diffex = c(-1, 1, 0)
                                        ),
      grouping_col = "gene",
      diffex_col   = "diffex"
    ),
    info = "Grouping column is absent - but default colname as arg"
  )

  expect_error(
    object = .get_consistent_diffexes(
      dframe       = tibble::data_frame(gene = letters[1:3],
                                        diffex = c(-1, 1, 0)
                                        ),
      grouping_col = "gene",
      diffex_col   = "NOT PRESENT"
    ),
    info = "Diffex column is absent"
  )

  expect_error(
    object = .get_consistent_diffexes(
      dframe       = tibble::data_frame(gene = letters[1:3],
                                        is.diffexed = c(-1, 1, 0)
                                        ),
      grouping_col = "gene",
      diffex_col   = "diffex"
    ),
    info = "Diffex column is absent - but default colname as arg"
  )

  expect_equal(
    object = .get_consistent_diffexes(
      dframe = tibble::data_frame(
        gene   = c("a", "b", "c"),
        diffex = c(1, 0, -1)
      ),
      grouping_col = "gene",
      diffex_col   = "diffex"
    ),
    expected = c("a", "c"),
    info = "Single gene tests"
  )

  expect_equal(
    object = .get_consistent_diffexes(
      dframe = tibble::data_frame(
        gene   = c("a", "a", "b", "c"),
        diffex = c(1, 1, 0, -1)
      ),
      grouping_col = "gene",
      diffex_col   = "diffex"
    ),
    expected = c("a", "c"),
    info = "One repeated gene"
  )

  expect_equal(
    object = .get_consistent_diffexes(
      dframe = tibble::data_frame(
        gene   = c("a", "a", "b", "b", "c", "c"),
        diffex = c(1,   0,   1,   -1,   0,  0)
      ),
      grouping_col = "gene",
      diffex_col   = "diffex"
    ),
    expected = character(0),
    info = "Repeated gene - but none are reproducibly diffexed"
  )

  expect_equal(
    object = .get_consistent_diffexes(
      dframe = tibble::data_frame(
        gene.id = c("a", "a", "b", "b", "c", "c"),
        diffex  = c(1,   1,   1,   -1,   1,  0)
      ),
      grouping_col = "gene.id",
      diffex_col   = "diffex"
    ),
    expected = "a",
    info = "Nonstandard gene column"
  )

  expect_equal(
    object = .get_consistent_diffexes(
      dframe = tibble::data_frame(
        gene      = c("a", "a", "b", "b", "c", "c"),
        diffexed  = c(1,   1,   1,   -1,   1,  0)
      ),
      grouping_col = "gene",
      diffex_col   = "diffexed"
    ),
    expected = "a",
    info = "Nonstandard diffex column"
  )

})
