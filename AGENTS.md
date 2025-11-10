# Repository Guidelines

## Project Structure & Module Organization
Core R APIs live in `R/` (e.g., `desurv_fit.R`, `cv_helpers.R`) and are exported through `NAMESPACE`. Native routines that back the hot paths are in `src/` (`functions.cpp`, associated objects, and platform `Makevars`). User-facing docs and pkg metadata live in `DESCRIPTION`, `man/`, and `README.Rmd`. Tests reside in `tests/testthat`, mirroring the function under test (`test-input_validation.R` today). Keep lightweight data fixtures under `tests/testthat/fixtures` if new cases are needed. Avoid editing `RcppExports.*` by hand—regenerate via `Rcpp::compileAttributes()`.

## Build, Test, and Development Commands
Use `devtools::load_all()` to iterate without reinstalling. Run `devtools::document()` before commits that touch roxygen comments so `man/` and `NAMESPACE` stay in sync. Package builds should pass `R CMD build .` followed by `R CMD check DeSurv_*.tar.gz` or the shortcut `devtools::check()` (which wraps build, test, lint). Native code can be recompiled quickly with `Rcpp::sourceCpp('src/functions.cpp')` when working outside the package context.

## Coding Style & Naming Conventions
Follow tidyverse-style R conventions: 2-space indent, snake_case objects, S3 methods as `generic.class`. Keep argument names descriptive and validate inputs via helpers in `validate_input.R`. For C++, prefer `camelCase` functions, `const auto&` where possible, and guard exported routines with thorough `// [[Rcpp::export]]` docs. Run `styler::style_pkg()` before large refactors to keep diffs clean.

## Testing Guidelines
Place new tests in `tests/testthat/test-<topic>.R`, grouping related expectations with `test_that()`. Mirror real data shapes by reusing synthetic fixtures and check both numeric accuracy and error messaging. Run `devtools::test()` locally before pushing; full CI parity requires `devtools::check()` because it re-runs the tests under multiple environments. Use `covr::report()` to confirm coverage stays at or above 85% for any PR touching statistical cores.

## Commit & Pull Request Guidelines
Commits should be small, written in the imperative, and scoped (e.g., "add cv helper for folds"). Reference issues with `Fixes #123` in the body when relevant. PRs need a summary of motivation, bullet list of changes, testing evidence (`devtools::check`, targeted scripts), and any screenshots/plots when visual behavior shifts. Tag reviewers familiar with the touched modules and highlight breaking changes early to keep release notes accurate.
