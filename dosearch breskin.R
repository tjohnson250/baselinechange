library(dosearch)

# Define the graph structure with bidirectional edge
# A = vaccine, W = pain, Y = disease, S = dropout
# W <-> Y represents unobserved confounding (U)
graph <- "
  A -> W
  A -> Y
  W <-> Y
  W -> S
"

# Define the query - we want the causal effect P(Y(a))
query <- "P(Y | do(A))"

# Define the data - what we can observe
# We observe A, W, Y when S=0 (not dropped out)
data <- "P(W | do(A))
  P(Y, W| do(A), S)"

# Run dosearch
result <- dosearch(
  data = data,
  query = query,
  graph = graph,
  selection_bias = "S"  # S is the selection variable
)

print(result)

# Expected from Breskin, at al.
# E(Y|do(a)) = Σ_w E(Y|W=w, S=0, A=a) × P(W=w|A=a)