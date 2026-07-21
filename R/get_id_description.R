#' Get Description and Ontology for GO/KEGG IDs
#'
#' Look up human-readable descriptions for Gene Ontology (GO) term IDs or
#' KEGG pathway IDs. The input is expected to be a single source at a time —
#' all GO IDs or all KEGG IDs, never mixed. GO lookups read the local
#' \code{GO.db} database (no network); KEGG lookups query \code{KEGGREST}
#' (requires an internet connection). The organism/reference source for KEGG
#' is inferred from the ID prefix, e.g. \code{hsa04110} queries
#' \emph{Homo sapiens}, \code{ath00010} \emph{Arabidopsis thaliana}, and
#' \code{map00010} the global reference pathway set.
#'
#' The returned columns follow the annotation-database schema used by
#' \code{\link{enrich_go}} and \code{\link{enrich_kegg}} (without the
#' \code{gene} column), so the result can be joined into a \code{go_db} or
#' \code{kegg_db} directly.
#'
#' @param ids Character vector of GO IDs (e.g. \code{"GO:0008150"}) or KEGG
#'   pathway IDs (e.g. \code{"hsa04110"}, \code{"map00010"}). Set to
#'   \code{NULL} to return the full mapping table for the source named by
#'   \code{type}.
#' @param type Character. How to interpret \code{ids}:
#'   \describe{
#'     \item{\code{"auto"} (default)}{Detect the source by format: if any ID
#'       matches \code{GO:#######} the vector is treated as GO, otherwise
#'       KEGG.}
#'     \item{\code{"go"}}{Treat the IDs as GO terms.}
#'     \item{\code{"kegg"}}{Treat the IDs as KEGG pathways.}
#'   }
#'   Must be \code{"go"} or \code{"kegg"} (not \code{"auto"}) when
#'   \code{ids = NULL}.
#' @param species KEGG organism code (e.g. \code{"hsa"}, \code{"ath"}) that
#'   overrides prefix inference. Only consulted for KEGG. If \code{NULL}
#'   (default), the organism is inferred from the shared prefix of the KEGG
#'   IDs; the reference prefix \code{map} selects the global pathway list.
#'   Ignored for GO.
#'
#' @return A \link[tibble]{tibble} whose columns depend on the source:
#'   \itemize{
#'     \item GO: \code{go_id}, \code{go_term}, \code{go_ontology} (BP/CC/MF)
#'     \item KEGG: \code{kegg_id}, \code{kegg_term}, \code{kegg_category}
#'       (\code{kegg_category} is always \code{NA} — KEGGREST does not
#'       provide pathway categories)
#'   }
#'   Input order is preserved in lookup mode. IDs with no match are kept with
#'   \code{NA} terms and reported via a warning. In full-table mode
#'   (\code{ids = NULL}) every record from the source is returned.
#'
#' @details KEGG pathway descriptions returned by KEGGREST carry an organism
#'   suffix for organism-specific queries (e.g.
#'   \code{"Cell cycle - Homo sapiens (human)"}); reference (\code{map})
#'   descriptions do not. The raw KEGG string is returned unchanged.
#'
#' @seealso \code{\link{enrich_go}}, \code{\link{enrich_kegg}}
#' @author Xiang LI <lixiang117423@gmail.com>
#' @export
#'
#' @examples
#' \dontrun{
#' # GO terms — local GO.db, no network needed
#' get_id_description(c("GO:0008150", "GO:0005737"))
#'
#' # KEGG pathways, organism inferred from the "hsa" prefix (network)
#' get_id_description(c("hsa04110", "hsa04151"))
#'
#' # KEGG reference pathways
#' get_id_description(c("map00010", "map00020"))
#'
#' # Force a source explicitly
#' get_id_description(c("0008150", "0005737"), type = "go")
#'
#' # Full mapping table
#' get_id_description(NULL, type = "go")
#' get_id_description(NULL, type = "kegg", species = "ath")
#' }
get_id_description <- function(ids = NULL, type = c("auto", "go", "kegg"),
                               species = NULL) {
  type <- match.arg(type)

  # Full-table mode: source must be explicit, nothing to infer from
  if (is.null(ids)) {
    if (type == "auto") {
      stop("When 'ids' is NULL, 'type' must be \"go\" or \"kegg\" (not \"auto\").")
    }
    if (type == "go") return(fetch_go_table())
    return(fetch_kegg_table(species))
  }

  if (!is.character(ids)) {
    stop("'ids' must be a character vector (or NULL for the full table).")
  }
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) {
    stop("'ids' is empty after dropping NA/empty strings.")
  }

  # Input is a single source (all GO or all KEGG); detect once.
  src <- if (type == "auto") {
    if (any(is_go_id(ids))) "go" else "kegg"
  } else {
    type
  }

  if (src == "go") build_go_result(ids) else build_kegg_result(ids, species)
}

# --- GO (local GO.db, offline) ----------------------------------------------

is_go_id <- function(x) {
  grepl("^GO:?[0-9]{7}$", x, ignore.case = TRUE)
}

normalize_go_ids <- function(x) {
  x <- toupper(trimws(x))
  sub("^(GO:?)?([0-9]{7})$", "GO:\\2", x)
}

# Generic internal frame: id, description, ontology
lookup_go <- function(go_ids) {
  # select() returns NA rows for individual invalid keys, but errors when
  # *none* of the keys are valid — fall back to NA rows so the caller's
  # unmatched-id warning can report them.
  res <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      GO.db::GO.db,
      keys = go_ids,
      columns = c("TERM", "ONTOLOGY"),
      keytype = "GOID"
    )),
    error = function(e) NULL
  )
  if (is.null(res)) {
    return(tibble::tibble(
      id = go_ids,
      description = NA_character_,
      ontology = NA_character_
    ))
  }
  tibble::tibble(id = res$GOID, description = res$TERM, ontology = res$ONTOLOGY)
}

build_go_result <- function(ids) {
  db <- lookup_go(unique(normalize_go_ids(ids)))
  out <- tibble::tibble(
    go_id = normalize_go_ids(ids),
    go_term = NA_character_,
    go_ontology = NA_character_
  )
  matched <- match(out$go_id, db$id)
  hit <- !is.na(matched)
  out$go_term[hit] <- db$description[matched[hit]]
  out$go_ontology[hit] <- db$ontology[matched[hit]]
  warn_unmatched(out$go_id, out$go_term)
  out
}

fetch_go_table <- function() {
  db <- lookup_go(AnnotationDbi::keys(GO.db::GO.db, "GOID"))
  tibble::tibble(go_id = db$id, go_term = db$description, go_ontology = db$ontology)
}

# --- KEGG (KEGGREST, network) -----------------------------------------------

kegg_clean_id <- function(x) {
  sub("^path:", "", trimws(x))
}

resolve_kegg_organism <- function(kegg_ids, species) {
  if (!is.null(species)) return(species)
  prefixes <- unique(sub("[0-9].*$", "", kegg_clean_id(kegg_ids)))
  if (length(prefixes) > 1) {
    stop("KEGG IDs have mixed prefixes (",
         paste(prefixes, collapse = ", "),
         "). Set 'species' or query one organism at a time.",
         call. = FALSE)
  }
  if (prefixes == "map") return(NA_character_)  # global reference pathways
  prefixes
}

kegg_pathway_list <- function(species) {
  tryCatch(
    if (is.null(species)) KEGGREST::keggList("pathway") else KEGGREST::keggList("pathway", species),
    error = function(e) {
      stop("KEGGREST query failed: ", e$message,
           ". Check the organism code / network connection.", call. = FALSE)
    }
  )
}

# Generic internal frame: id, description
lookup_kegg <- function(kegg_ids, species) {
  org <- resolve_kegg_organism(kegg_ids, species)
  raw <- if (is.na(org)) kegg_pathway_list(NULL) else kegg_pathway_list(org)
  tbl_id <- kegg_clean_id(names(raw))
  matched <- match(kegg_clean_id(kegg_ids), tbl_id)
  tibble::tibble(
    id = kegg_ids,
    description = ifelse(is.na(matched), NA_character_, unname(raw)[matched])
  )
}

build_kegg_result <- function(ids, species) {
  db <- lookup_kegg(unique(ids), species)
  out <- tibble::tibble(
    kegg_id = ids,
    kegg_term = NA_character_,
    kegg_category = NA_character_  # KEGGREST keggList does not provide categories
  )
  matched <- match(out$kegg_id, db$id)
  hit <- !is.na(matched)
  out$kegg_term[hit] <- db$description[matched[hit]]
  warn_unmatched(out$kegg_id, out$kegg_term)
  out
}

fetch_kegg_table <- function(species) {
  raw <- kegg_pathway_list(species)
  tibble::tibble(
    kegg_id = kegg_clean_id(names(raw)),
    kegg_term = unname(raw),
    kegg_category = NA_character_
  )
}

# --- shared -----------------------------------------------------------------

warn_unmatched <- function(ids, terms) {
  missing_ids <- unique(ids[is.na(terms)])
  if (length(missing_ids)) {
    preview <- paste(head(missing_ids, 10), collapse = ", ")
    warning("No description found for ", length(missing_ids), " id(s): ",
            preview, if (length(missing_ids) > 10) " ...", call. = FALSE)
  }
  invisible(NULL)
}
