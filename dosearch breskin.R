library(dosearch)

# Deriving the g-formula in Breskin et al.: 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5882566/#R4
# Breskin A, Cole SR, Hudgens MG. A practical example demonstrating the utility
# of single-world intervention graphs. Epidemiology. 2018 May 1;29(3):e20-1.

# Define the graph structure with bidirectional edge
# A = treatment, Y = outcome, S = dropout
# W <-> Y represents unobserved confounding (U)
graph <- "
  A -> W
  A -> Y
  W <-> Y
  W -> S
"

# Define the query - we want the causal effect P(Y(a))
query <- "P(Y | do(A))"

# Define the data distributions - what we can observe
# We observe W given do(A) and
# Y and W given do(A) and S (defined in the dosearch call as the selection variable)
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
#> \sum_{W}\left(p(W|do(A))p(Y|do(A),W,S)\right)  

# Expected from Breskin, at al.
# E(Y|do(a)) = Σ_w E(Y|W=w, S=0, A=a) × P(W=w|A=a)