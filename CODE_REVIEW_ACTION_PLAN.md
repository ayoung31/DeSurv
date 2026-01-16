# DeSurv Code Review: Validated Findings & Action Plan

**Date:** 2026-01-07
**Branch:** `tie_handling`
**Status:** Ready for implementation

---

## Executive Summary

After initial code review and validation pass, **9 confirmed issues** were identified in the DeSurv package. This document consolidates the validated findings and provides a prioritized action plan.

| Severity | Count |
|----------|-------|
| Important | 5 |
| Minor | 4 |

---

## Part 1: Confirmed Issues

### Important Issues

#### 1. Missing Nonnegativity Validation for X
**Files:** `R/validate_input.R`
**Impact:** High

X is not checked for negative values. Negative input violates NMF assumptions and can break initialization (e.g., `runif(n, 0, max(X))` with a negative max in `init()`), or cause invalid multiplicative updates.

**Current behavior:** Allows any numeric values including negatives.

**Expected behavior:** Reject X with any negative values since NMF requires non-negative input.

---

#### 2. CV Mean/SE Calculation Inconsistency
**Files:** `R/cv_warm.R`, `R/cv_cold.R`
**Impact:** High

`mean_cindex` returns `NA` if **any** fold is non-finite, while `se_cindex` drops non-finite values first. This produces logically inconsistent summaries where `mean_cindex = NA` but `se_cindex` has a finite value.

**Current behavior:**
```r
# Mean - returns NA if any non-finite
FUN = function(x) {
  if (anyNA(x) || any(!is.finite(x))) return(NA_real_)
  mean(x)
}

# SE - filters out non-finite first
se_fun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1L) return(NA_real_)
  stats::sd(x) / sqrt(length(x))
}
```

**Expected behavior:** Both functions should handle non-finite values consistently (either both filter or both return NA).

---

#### 3. predict() Does Not Validate Finite Values in newdata
**Files:** `R/predict_methods.R`
**Impact:** Medium-High

`newdata` can contain `NA/NaN/Inf` and is not rejected, which propagates non-finite values into predictions silently.

**Current behavior:** Only checks `is.numeric()`, allows non-finite values.

**Expected behavior:** Reject or warn when newdata contains non-finite values.

---

#### 4. H Output Order Mismatch Due to Time Sorting (NEW)
**Files:** `src/functions.cpp`
**Impact:** High

`optimize_loss_cpp()` sorts samples by time (line 520-524) and reorders H accordingly, but **never restores the original order** before returning. The returned H columns no longer align with the input sample order.

**Current behavior:**
```cpp
// Sort by increasing time once
arma::uvec tOrder = arma::sort_index(y);
y = y.elem(tOrder);
d = d.elem(tOrder);
X = X.cols(tOrder);
H = H.cols(tOrder);  // H is reordered...
// ... optimization happens ...
// H is returned in time-sorted order, not original order
```

**Expected behavior:** Either restore H to original column order before returning, or clearly document that H columns are in time-sorted order.

---

#### 5. All-Zero X Yields Division by Zero (NEW)
**Files:** `src/functions.cpp`, `R/validate_input.R`
**Impact:** Medium-High

`Xnorm = sum(X²)` can be 0 (e.g., all-zero data), which causes divisions by zero in `calc_loss_cpp()` and update steps.

**Current behavior:**
```cpp
double Xnorm = arma::accu(arma::square(X));  // Can be 0
// Later used in divisions:
double nmf_loss = arma::accu(arma::square(X - W * H))/(2.0 * Xnorm);  // Division by 0!
```

**Expected behavior:** Validate that `sum(X²) > 0` or handle the edge case gracefully.

---

### Minor Issues

#### 6. Alpha Clipping Contradicts Documentation
**Files:** `R/validate_input.R`, `R/desurv_fit.R`
**Impact:** Low

Alpha values >= 1 are clipped to `1 - .Machine$double.eps`, while documentation says alpha is in [0, 1] inclusive.

**Fix:** Update documentation to state `[0, 1)` or change validation to accept alpha=1.

---

#### 7. `max_iter_beta` Is Undocumented
**Files:** `R/desurv_fit.R`
**Impact:** Low

Parameter is exposed in the function signature (line 211) but missing from roxygen documentation.

**Fix:** Add `@param max_iter_beta` to roxygen comments.

---

#### 8. init() Hides Failure Context in Parallel Mode
**Files:** `R/init.R`
**Impact:** Low

Parallel path suppresses per-init warnings (`verbose_inner = FALSE`) and the final error only states that all initializations failed, making debugging difficult.

**Fix:** Collect and report error messages from failed initializations.

---

#### 9. BO Default Exploration Is Disabled
**Files:** `R/desurv_bayesopt.R`, `R/desurv_cv_bayesopt.R`
**Impact:** Low

`exploration_weight` defaults to 0, so Expected Improvement is fully exploitative unless the user overrides it. Standard default is typically 0.01.

**Fix:** Consider changing default to 0.01 or document the purely exploitative behavior.

---

## Part 2: Test Coverage Gaps

The following tests should be added to cover the confirmed issues:

| Test | File | Issue Covered |
|------|------|---------------|
| Reject negative X values | `test-desurv_data.R` | #1 |
| Non-finite newdata in predict() | `test-predict-preprocess.R` | #3 |
| All-zero X input (Xnorm=0) | `test-desurv_data.R` | #5 |
| H column order matches input after fit | `test-fit-workflow.R` | #4 |
| alpha=1 clipping behavior | `test-desurv_data.R` | #6 |
| CV mean/SE consistency with NA folds | `test-cv-folds.R` | #2 |

---

## Part 3: Action Plan

### Phase 1: Critical Input Validation (Priority: Immediate)

| Task | File | Effort |
|------|------|--------|
| Add `any(X < 0)` check in `.validate_desurv_data()` | `R/validate_input.R` | Small |
| Add `sum(X²) > 0` check or handle Xnorm=0 | `R/validate_input.R` | Small |
| Add finite value check in `.desurv_prepare_newdata()` | `R/predict_methods.R` | Small |
| Add tests for all three validations | `tests/testthat/` | Small |

**Estimated effort:** 1-2 hours

---

### Phase 2: Fix CV Inconsistency (Priority: High)

| Task | File | Effort |
|------|------|--------|
| Align mean calculation with SE calculation (filter non-finite first) | `R/cv_warm.R` | Small |
| Same fix in cold engine | `R/cv_cold.R` | Small |
| Add test for consistent mean/SE with NA folds | `tests/testthat/test-cv-folds.R` | Small |

**Estimated effort:** 30 minutes

---

### Phase 3: Fix H Column Order (Priority: High)

| Task | File | Effort |
|------|------|--------|
| Option A: Restore H to original order before returning from `optimize_loss_cpp()` | `src/functions.cpp` | Medium |
| Option B: Document that H is in time-sorted order | `R/desurv_fit.R` | Small |
| Add test verifying H column alignment | `tests/testthat/test-fit-workflow.R` | Small |

**Recommended:** Option A (restore order) for API consistency.

**Estimated effort:** 1 hour

---

### Phase 4: Documentation & Minor Fixes (Priority: Low)

| Task | File | Effort |
|------|------|--------|
| Add `@param max_iter_beta` documentation | `R/desurv_fit.R` | Trivial |
| Update alpha documentation to say `[0, 1)` | `R/desurv_fit.R`, `man/` | Trivial |
| Consider changing BO `exploration_weight` default to 0.01 | `R/desurv_bayesopt.R` | Trivial |
| Improve init() error messages in parallel mode | `R/init.R` | Small |
| Run `devtools::document()` to regenerate docs | - | Trivial |

**Estimated effort:** 30 minutes

---

## Part 4: Implementation Checklist

```
Phase 1: Input Validation
[ ] Add X >= 0 check in validate_input.R
[ ] Add Xnorm > 0 check in validate_input.R
[ ] Add finite check in predict_methods.R
[ ] Add tests for negative X
[ ] Add tests for all-zero X
[ ] Add tests for non-finite newdata

Phase 2: CV Consistency
[ ] Fix mean calculation in cv_warm.R
[ ] Fix mean calculation in cv_cold.R
[ ] Add test for mean/SE consistency

Phase 3: H Column Order
[ ] Restore H order in functions.cpp OR document behavior
[ ] Add test for H column alignment

Phase 4: Documentation
[ ] Document max_iter_beta
[ ] Update alpha range in docs
[ ] Review exploration_weight default
[ ] Improve init() error messages
[ ] Run devtools::document()
[ ] Run devtools::check()
```

---

## Part 5: Verification

After implementing fixes, verify:

```r
# All tests pass
devtools::test()

# Package check passes
devtools::check()

# Manual verification of key fixes
X_neg <- matrix(c(-1, 2, 3, 4), nrow = 2)
desurv_data(X_neg, y, d, k = 2)  # Should error

X_zero <- matrix(0, nrow = 5, ncol = 10)
desurv_data(X_zero, y, d, k = 2)  # Should error or warn

fit <- desurv_fit(X, y, d, k = 2, ...)
all(colnames(fit$H) == colnames(X))
# OR H columns match original sample order
```

---

## Summary

| Phase | Issues Addressed | Effort | Priority |
|-------|------------------|--------|----------|
| 1 | #1, #3, #5 | 1-2 hours | Immediate |
| 2 | #2 | 30 min | High |
| 3 | #4 | 1 hour | High |
| 4 | #6, #7, #8, #9 | 30 min | Low |

**Total estimated effort:** 3-4 hours

The most impactful fixes are the input validation (Phase 1) and H column order (Phase 3), which prevent silent data corruption. The CV consistency fix (Phase 2) improves reliability of model selection.

---

*Generated from validated code review findings.*
