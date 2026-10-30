# Contributing to bioRtools

Contributions are welcome through GitHub issues and pull requests. By
participating, you agree to follow the project [Code of Conduct](CODE_OF_CONDUCT.md).

## Report a bug

Open an issue at <https://github.com/lixiang117423/bioRtools/issues> and include:

- the installed `bioRtools` version;
- the output of `sessionInfo()`;
- a minimal reproducible example using public or synthetic data;
- the command that was run, the expected result, and the observed result.

Do not include confidential research data, credentials, or personal information
in a public issue.

## Request a feature

Describe the research problem, the input and output you expect, and which
existing function or workflow is closest to the request. Proposed additions
should provide reusable research value rather than duplicate an existing
function with a different name.

## Development workflow

1. Fork the repository and create a branch from `main`.
2. Install development dependencies with `pak::pak("lixiang117423/bioRtools")`.
3. Add or update tests in `tests/testthat/`.
4. Update roxygen documentation and run `devtools::document()`.
5. Run `devtools::test()` and `devtools::check()`.
6. Open a pull request describing the behaviour change and the tests performed.

For a focused local check:

```r
devtools::document()
devtools::test()
devtools::check(args = c("--no-manual"), error_on = "warning")
```

## Coding conventions

- Use snake_case for function names, arguments, and output columns.
- Use explicit package namespaces such as `dplyr::filter()`.
- Return objects rather than relying on printing or plotting as a side effect.
- Validate sample identifiers and alignment whenever two tables are combined.
- Use English error messages that identify the invalid argument and expected
  condition.
- Add a test for every bug fix and every exported function or changed result
  contract.

## Pull-request checklist

- [ ] The change has a focused research or maintenance purpose.
- [ ] Documentation and examples are updated.
- [ ] Tests cover the new or corrected behaviour.
- [ ] `devtools::test()` passes.
- [ ] `R CMD check` has no errors or warnings.
- [ ] NEWS or `CHANGELOG.md` is updated for user-visible changes.
- [ ] No confidential data, credentials, or absolute local paths are committed.

## Support

Usage questions belong in GitHub Discussions when enabled, or in an issue with
the `question` label. Security-sensitive reports should follow
[SUPPORT.md](SUPPORT.md).
