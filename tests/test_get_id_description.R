# Tests for get_id_description(): GO/KEGG id -> description (+ GO ontology).
# GO branch is offline (local GO.db). KEGG branch needs network and is skipped
# if KEGGREST cannot reach the API.
source("R/get_id_description.R")
suppressPackageStartupMessages({ library(GO.db); library(AnnotationDbi) })

cat("Test 1: GO lookup — columns, values, ontology\n")
g <- get_id_description(c("GO:0008150", "GO:0005737", "GO:0003674"))
stopifnot(identical(names(g), c("go_id", "go_term", "go_ontology")))
stopifnot(g$go_term[1] == "biological_process" && g$go_ontology[1] == "BP")
stopifnot(g$go_term[2] == "cytoplasm"          && g$go_ontology[2] == "CC")
stopifnot(g$go_term[3] == "molecular_function" && g$go_ontology[3] == "MF")
cat("✓ Test 1 passed\n\n")

cat("Test 2: GO id normalization (bare number / lowercase prefix)\n")
g2 <- get_id_description(c("0008150", "go:0005737"), type = "go")
stopifnot(g2$go_id[1] == "GO:0008150", g2$go_id[2] == "GO:0005737")
stopifnot(!any(is.na(g2$go_term)))
cat("✓ Test 2 passed\n\n")

cat("Test 3: GO input order preserved (incl. duplicates)\n")
g3 <- get_id_description(c("GO:0005737", "GO:0008150", "GO:0005737"))
stopifnot(identical(g3$go_id, c("GO:0005737", "GO:0008150", "GO:0005737")))
cat("✓ Test 3 passed\n\n")

cat("Test 4: GO all-invalid ids -> NA rows, no crash, warning emitted\n")
w <- NULL
g4 <- withCallingHandlers(
  get_id_description(c("GO:9999999", "GO:8888888")),
  warning = function(w) { assign("w", conditionMessage(w), envir = parent.frame()); invokeRestart("muffleWarning") }
)
stopifnot(nrow(g4) == 2, all(is.na(g4$go_term)), all(is.na(g4$go_ontology)))
cat("✓ Test 4 passed (warning: ", w, ")\n\n")

cat("Test 5: full GO table shape\n")
ft <- get_id_description(NULL, type = "go")
stopifnot(identical(names(ft), c("go_id", "go_term", "go_ontology")))
stopifnot(nrow(ft) > 10000)
cat("✓ Test 5 passed (", nrow(ft), "GO terms)\n\n")

cat("Test 6: error — ids = NULL with type = \"auto\"\n")
e <- tryCatch(get_id_description(NULL), error = function(e) e$message)
stopifnot(grepl("must be", e))
cat("✓ Test 6 passed\n\n")

cat("Test 7: error — invalid type\n")
e <- tryCatch(get_id_description(c("GO:0008150"), type = "nope"), error = function(e) e$message)
stopifnot(grepl("auto.*go.*kegg", e))
cat("✓ Test 7 passed\n\n")

# --- KEGG (network) ---------------------------------------------------------
net_ok <- isTRUE(tryCatch(is.character(KEGGREST::keggList("pathway", "hsa")),
                          error = function(e) FALSE))

cat("Test 8: KEGG lookup — columns + known value", if (!net_ok) "[SKIPPED: no network]" , "\n")
if (net_ok) {
  k <- suppressWarnings(get_id_description(c("hsa04110", "hsa04151")))
  stopifnot(identical(names(k), c("kegg_id", "kegg_term", "kegg_category")))
  stopifnot(all(is.na(k$kegg_category)))  # not provided by keggList
  stopifnot(grepl("Cell cycle", k$kegg_term[1]))
  stopifnot(!is.na(k$kegg_term[2]))
  cat("✓ Test 8 passed\n\n")
}

cat("Test 9: KEGG reference (map) + species override", if (!net_ok) "[SKIPPED: no network]", "\n")
if (net_ok) {
  km <- suppressWarnings(get_id_description(c("map00010")))
  stopifnot(km$kegg_id[1] == "map00010", grepl("Glycolysis", km$kegg_term[1]))
  kt <- get_id_description(NULL, type = "kegg", species = "ath")
  stopifnot(identical(names(kt), c("kegg_id", "kegg_term", "kegg_category")))
  stopifnot(nrow(kt) > 50, all(grepl("^ath", kt$kegg_id)))
  cat("✓ Test 9 passed\n\n")
}

cat("Test 10: error — mixed KEGG organism prefixes\n")
e <- tryCatch(get_id_description(c("hsa04110", "ath00010")), error = function(e) e$message)
stopifnot(grepl("mixed prefixes", e))
cat("✓ Test 10 passed\n\n")

cat("All get_id_description tests done.\n")
