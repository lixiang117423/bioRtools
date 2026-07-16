# Tests for predict_func() (FAPROTAX) + match_words().
source("R/parse_faprotax.R")
source("R/predict_func.R")

cat("Test 1: parse_faprotax — ~90 functions, methanotrophy clades\n")
parsed <- parse_faprotax()
stopifnot(length(parsed$func_tax) >= 85 && length(parsed$func_tax) <= 105)
stopifnot("methanotrophy" %in% names(parsed$func_tax))
stopifnot(any(grepl("Methylococcaceae", parsed$func_tax$methanotrophy)))
stopifnot(!is.na(parsed$ver))
cat("✓ Test 1 passed (", length(parsed$func_tax), "functions, ver", parsed$ver, ")\n\n")

cat("Test 2: match_words boundary semantics (port of collapse_table.py)\n")
# matcher is case-sensitive; caller normalizes case. Use lowercase both sides.
stopifnot( match_words("*methylococcaceae*", "bacteria;proteobacteria;methylococcaceae"))
stopifnot(!match_words("*europaea*", "nitroeuropaea"))                       # word boundary
stopifnot( match_words("*proteobacteria*europaea*", "bacteria;proteobacteria;nitrosomonas;europaea"))
stopifnot(!match_words("*proteobacteria*europaea*", "proteobacteria;nitroeuropaea"))  # europaea not a word here
stopifnot( match_words("*e*coli*", "e;coli"))                                # two words
stopifnot(!match_words("*e*coli*", "escherichia;coli"))                      # 'e' glued to 'scherichia'
cat("✓ Test 2 passed\n\n")

cat("Test 3: predict_func assigns Methylococcaceae -> methanotrophy\n")
taxonomy <- data.frame(
  Kingdom = c("Bacteria","Bacteria","Archaea"),
  Phylum  = c("Proteobacteria","Proteobacteria","Euryarchaeota"),
  Class   = c("Gammaproteobacteria","Betaproteobacteria","Methanomicrobia"),
  Order   = c("Methylococcales","Nitrosomonadales","Methanosarcinales"),
  Family  = c("Methylococcaceae","Nitrosomonadaceae","Methanosarcinaceae"),
  Genus   = c("Methylosinus","Nitrosomonas","Methanosarcina"),
  Species = c("Methylosinus_trichosporium","Nitrosomonas_europaea","Methanosarcina_mazei"),
  row.names = c("ASV1","ASV2","ASV3"), stringsAsFactors = FALSE)
otu <- matrix(c(100,50,20, 80,40,10, 60,30,5), nrow = 3, ncol = 3,
              dimnames = list(c("ASV1","ASV2","ASV3"), c("S1","S2","S3")))
res <- predict_func(otu, taxonomy, verbose = FALSE)
stopifnot(res$asv_func["ASV1","methanotrophy"] == 1L)
stopifnot(all(res$asv_func["ASV2","methanotrophy"] == 0L,
              res$asv_func["ASV3","methanotrophy"] == 0L))
cat("✓ Test 3 passed\n\n")

cat("Test 4: func_table == t(t(m) %*% asv_func); methanotrophy row = ASV1 counts\n")
stopifnot(identical(dim(res$func_table), c(length(colSums(res$asv_func)), 3L)))
stopifnot(isTRUE(all.equal(as.numeric(res$func_table["methanotrophy", ]),
                           as.numeric(otu["ASV1", ]))))
# global consistency
recomputed <- t(t(otu) %*% res$asv_func)
stopifnot(isTRUE(all.equal(res$func_table, recomputed)))
cat("✓ Test 4 passed\n\n")

cat("Test 5: errors\n")
stopifnot(inherits(try(predict_func(otu, taxonomy, db = "nope"), silent = TRUE), "try-error"))
bad_tax <- taxonomy; rownames(bad_tax) <- c("X1","X2","X3")
stopifnot(inherits(try(predict_func(otu, bad_tax, verbose = FALSE), silent = TRUE), "try-error"))
cat("✓ Test 5 passed\n\n")

cat("=====================================\n")
cat("All predict_func tests passed! ✓\n")
cat("=====================================\n")
