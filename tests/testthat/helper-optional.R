can_load_namespace <- function(package) {
  isTRUE(tryCatch(
    {
      loadNamespace(package)
      TRUE
    },
    error = function(e) FALSE
  ))
}

can_load_namespaces <- function(packages) {
  all(vapply(packages, can_load_namespace, logical(1)))
}

skip_if_namespace_unavailable <- function(package) {
  if (!can_load_namespace(package)) {
    testthat::skip(sprintf("{%s} is not available", package))
  }
}
