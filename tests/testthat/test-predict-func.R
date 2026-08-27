# Migrated from tests/test_predict_func.R (bare script).

match_words <- bioRtools:::match_words
parse_faprotax <- bioRtools:::parse_faprotax

make_taxonomy <- function() {
  taxonomy <- data.frame(
    Kingdom = c("Bacteria", "Bacteria", "Archaea"),
    Phylum  = c("Proteobacteria", "Proteobacteria", "Euryarchaeota"),
    Class   = c("Gammaproteobacteria", "Betaproteobacteria", "Methanomicrobia"),
    Order   = c("Methylococcales", "Nitrosomonadales", "Methanosarcinales"),
    Family  = c("Methylococcaceae", "Nitrosomonadaceae", "Methanosarcinaceae"),
    Genus   = c("Methylosinus", "Nitrosomonas", "Methanosarcina"),
    Species = c("Methylosinus_trichosporium", "Nitrosomonas_europaea", "Methanosarcina_mazei"),
    row.names = c("ASV1", "ASV2", "ASV3"), stringsAsFactors = FALSE
  )
  taxonomy
}

make_otu <- function() {
  matrix(c(100, 50, 20, 80, 40, 10, 60, 30, 5), nrow = 3, ncol = 3,
         dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2", "S3")))
}

test_that("parse_faprotax — ~90 functions, methanotrophy clades", {
  parsed <- parse_faprotax()
  expect_true(length(parsed$func_tax) >= 85 && length(parsed$func_tax) <= 105)
  expect_true("methanotrophy" %in% names(parsed$func_tax))
  expect_true(any(grepl("Methylococcaceae", parsed$func_tax$methanotrophy)))
  expect_false(is.na(parsed$ver))
})

test_that("match_words boundary semantics (port of collapse_table.py)", {
  expect_true(match_words("*methylococcaceae*", "bacteria;proteobacteria;methylococcaceae"))
  expect_false(match_words("*europaea*", "nitroeuropaea"))
  expect_true(match_words("*proteobacteria*europaea*", "bacteria;proteobacteria;nitrosomonas;europaea"))
  expect_false(match_words("*proteobacteria*europaea*", "proteobacteria;nitroeuropaea"))
  expect_true(match_words("*e*coli*", "e;coli"))
  expect_false(match_words("*e*coli*", "escherichia;coli"))
})

test_that("predict_func assigns Methylococcaceae -> methanotrophy", {
  otu <- make_otu()
  taxonomy <- make_taxonomy()
  res <- predict_func(otu, taxonomy, verbose = FALSE)
  expect_equal(res$asv_func["ASV1", "methanotrophy"], 1L)
  expect_true(all(res$asv_func["ASV2", "methanotrophy"] == 0L,
                  res$asv_func["ASV3", "methanotrophy"] == 0L))
})

test_that("func_table == t(t(m) %*% asv_func); methanotrophy row = ASV1 counts", {
  otu <- make_otu()
  taxonomy <- make_taxonomy()
  res <- predict_func(otu, taxonomy, verbose = FALSE)
  expect_identical(dim(res$func_table), c(length(colSums(res$asv_func)), 3L))
  expect_equal(as.numeric(res$func_table["methanotrophy", ]),
               as.numeric(otu["ASV1", ]))
  recomputed <- t(t(otu) %*% res$asv_func)
  expect_equal(res$func_table, recomputed)
})

test_that("errors", {
  otu <- make_otu()
  taxonomy <- make_taxonomy()
  expect_error(predict_func(otu, taxonomy, db = "nope"))
  bad_tax <- taxonomy
  rownames(bad_tax) <- c("X1", "X2", "X3")
  expect_error(predict_func(otu, bad_tax, verbose = FALSE))
})
