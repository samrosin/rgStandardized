# The inverse logit function
inv.logit <- function(x){1/(1 + exp(-x))}

# truncate a number into the range [0, 1]
truncate_01 <- function(x) {min(1, max(0, x))}
