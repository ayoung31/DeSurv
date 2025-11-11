make_fixture_dataset <- function(
    p = 8,
    n = 30,
    k = 2,
    seed = 42
) {
  set.seed(seed)
  X <- matrix(rexp(p * n), nrow = p, ncol = n)
  rownames(X) <- paste0("gene", seq_len(p))
  y <- rexp(n) + 0.1
  d <- rbinom(n, 1, 0.6)

  # ensure at least a few events for every fold
  d[seq_len(min(5, n))] <- 1
  d[(seq_len(min(5, n)) + 5)] <- 0

  list(
    X = X,
    y = y,
    d = d,
    k = k,
    p = p,
    n = n
  )
}
