subsets <- function(n, r, v = 1:n, set = TRUE) {
################################################
        if(r < 0 || r > n) stop("invalid r for this n")
        if(set) {
                v <- unique(sort(v))
                if (length(v) < n) stop("too few different elements")
        }
        v0 <- vector(mode(v), 0)
        sub <- function(n, r, v) { ## Inner workhorse
                if(r == 0) v0 else
                if(r == n) matrix(v, 1, n) else
                rbind(cbind(    v[1],
                                Recall(n-1, r-1, v[-1])),
                                Recall(n-1, r, v[-1]))
        }
        sub(n, r, v[1:n])
}
