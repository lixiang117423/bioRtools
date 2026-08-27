# Migrated from tests/test_write_data_by_group.R (bare script).

make_grouped_df <- function() {
  df <- tibble::tibble(
    CHROM = c(rep("chr1", 3), rep("chr2", 2)),
    pos = c(100, 200, 300, 400, 500),
    p_value = c(0.1, 0.02, 0.5, 0.001, 0.9)
  )
  df
}

out_dir <- file.path(tempdir(), "write_data_by_group_test")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

test_that("single-column group, xlsx", {
  df <- make_grouped_df()
  df %>%
    dplyr::group_by(CHROM) %>%
    write_data_by_group(file.path(out_dir, "result_{CHROM}.xlsx"))
  f1 <- file.path(out_dir, "result_chr1.xlsx")
  f2 <- file.path(out_dir, "result_chr2.xlsx")
  expect_true(file.exists(f1))
  expect_true(file.exists(f2))
  expect_equal(nrow(read_data(f1)), 3)
  expect_equal(nrow(read_data(f2)), 2)
})

test_that("multi-column group, csv", {
  df <- make_grouped_df()
  df$TYPE <- c("A", "A", "B", "B", "B")
  df %>%
    dplyr::group_by(CHROM, TYPE) %>%
    write_data_by_group(file.path(out_dir, "{CHROM}_{TYPE}.csv"))
  expect_true(file.exists(file.path(out_dir, "chr1_A.csv")))
  expect_true(file.exists(file.path(out_dir, "chr2_B.csv")))
})

test_that("returns input invisibly", {
  df <- make_grouped_df()
  res <- df %>%
    dplyr::group_by(CHROM) %>%
    write_data_by_group(file.path(out_dir, "p_{CHROM}.xlsx"))
  expect_true(dplyr::is_grouped_df(res))
})

test_that("error on ungrouped data", {
  df <- make_grouped_df()
  expect_error(write_data_by_group(df, file.path(out_dir, "x_{CHROM}.xlsx")))
})
