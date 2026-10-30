if (!can_load_namespaces(c("AnnotationDbi", "GO.db", "KEGGREST"))) {
  test_that("GO/KEGG annotation integration tests are skipped", {
    testthat::skip("AnnotationDbi, GO.db, or KEGGREST is not available")
  })
} else {

# Migrated from tests/test_get_id_description.R (bare script).
# get_id_description(): GO/KEGG id -> description (+ GO ontology).
# GO branch is offline (local GO.db). KEGG branch needs network and is
# skipped if KEGGREST cannot reach the API.

net_ok <- isTRUE(tryCatch(
  is.character(KEGGREST::keggList("pathway", "hsa")),
  error = function(e) FALSE
))

test_that("GO lookup — columns, values, ontology", {
  g <- get_id_description(c("GO:0008150", "GO:0005737", "GO:0003674"))

  expect_identical(names(g), c("go_id", "go_term", "go_ontology"))
  expect_equal(g$go_term[1], "biological_process")
  expect_equal(g$go_ontology[1], "BP")
  expect_equal(g$go_term[2], "cytoplasm")
  expect_equal(g$go_ontology[2], "CC")
  expect_equal(g$go_term[3], "molecular_function")
  expect_equal(g$go_ontology[3], "MF")
})

test_that("GO id normalization (bare number / lowercase prefix)", {
  g2 <- get_id_description(c("0008150", "go:0005737"), type = "go")

  expect_equal(g2$go_id[1], "GO:0008150")
  expect_equal(g2$go_id[2], "GO:0005737")
  expect_false(any(is.na(g2$go_term)))
})

test_that("GO input order preserved (incl. duplicates)", {
  g3 <- get_id_description(c("GO:0005737", "GO:0008150", "GO:0005737"))

  expect_identical(g3$go_id, c("GO:0005737", "GO:0008150", "GO:0005737"))
})

test_that("GO all-invalid ids -> NA rows, no crash", {
  g4 <- withCallingHandlers(
    get_id_description(c("GO:9999999", "GO:8888888")),
    warning = function(w) invokeRestart("muffleWarning")
  )

  expect_equal(nrow(g4), 2)
  expect_true(all(is.na(g4$go_term)))
  expect_true(all(is.na(g4$go_ontology)))
})

test_that("full GO table shape", {
  ft <- get_id_description(NULL, type = "go")

  expect_identical(names(ft), c("go_id", "go_term", "go_ontology"))
  expect_gt(nrow(ft), 10000)
})

test_that("error — ids = NULL with type = 'auto'", {
  expect_error(get_id_description(NULL), "must be")
})

test_that("error — invalid type", {
  expect_error(get_id_description(c("GO:0008150"), type = "nope"),
               "auto.*go.*kegg")
})

test_that("KEGG lookup — columns + known value", {
  skip_if(!net_ok, "no network / KEGG API unreachable")

  k <- suppressWarnings(get_id_description(c("hsa04110", "hsa04151")))
  expect_identical(names(k), c("kegg_id", "kegg_term", "kegg_category"))
  expect_true(all(is.na(k$kegg_category)))   # not provided by keggList
  expect_match(k$kegg_term[1], "Cell cycle")
  expect_false(is.na(k$kegg_term[2]))
})

test_that("KEGG reference (map) + species override", {
  skip_if(!net_ok, "no network / KEGG API unreachable")

  km <- suppressWarnings(get_id_description(c("map00010")))
  expect_equal(km$kegg_id[1], "map00010")
  expect_match(km$kegg_term[1], "Glycolysis")

  kt <- get_id_description(NULL, type = "kegg", species = "ath")
  expect_identical(names(kt), c("kegg_id", "kegg_term", "kegg_category"))
  expect_gt(nrow(kt), 50)
  expect_true(all(grepl("^ath", kt$kegg_id)))
})

test_that("error — mixed KEGG organism prefixes", {
  expect_error(get_id_description(c("hsa04110", "ath00010")),
               "mixed prefixes")
})
}
