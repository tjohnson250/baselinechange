library(causaleffect)
library(igraph)

# Define DAG
g <- graph_from_literal(X -+ Y, Z -+ X, Z -+ Y, W -+ Z)
# Check if P(Y|do(X=x)) is identifiable
result <- causal.effect("Y", "X", G = g)
print(result)  #

# Define DAG
# g <- graph_from_literal(A -+ W, A -+ Y, W -+ S, U -+ W, U -+ Y)
## W and Y are connected by U, a latent variable, but causal.effect expects U to be replaced
## with an undirected line as follows:
g <- graph_from_literal(A -+ W, A -+ Y, W -+ S, W -+ Y, Y -+ W, simplify = FALSE)
# Note that U is unobserved by setting edges between W and Y
g <- set_edge_attr(graph = g, name = "description", index = c(4, 5), value = "U")
# A selected variable must have its vertex description set to S
# Use V(g) to find vertex indices, here S is 4
g <- set.vertex.attribute(graph = g, 
                          name = "description", index = 4, value = "S")
# Check if P(Y|do(A=a)) is identifiable
result <- causal.effect("Y", "A", G = g, simp = TRUE)
print(result)  #

result <- causal.effect("Y", "A", "S", G = g, simp = TRUE)
print(result)

result <- causal.effect(c("Y", "S"), "A", G = g, simp = TRUE)
print(result)


# Step 1: Build the graph exactly as in the paper
g <- graph.formula(A -+ W -+ S, A -+ Y, U -+ W, U -+ Y)
g <- set_edge_attr(graph = g, name = "description", index = c(4,5), value = "U")
                   #index = E(g)[.from("U")], value = "U")


# Step 2: Visualize to verify structure
plot(g, 
     vertex.color = ifelse(V(g)$name == "U", "red", "lightblue"),
     vertex.label.cex = 1.2,
     edge.curved = 0.1,
     main = "Vaccine Trial DAG")

# Step 3: Check without selection first
no_selection <- causal.effect("Y", "A", G = g)
cat("Without selection:\n")
print(no_selection)

# Step 4: Check with selection
with_selection <- causal.effect("Y", "A", "S", G = g)
cat("With selection:\n")
print(with_selection)

# Step 5: Verify the confounding structure
cat("Is there a path from U to Y?\n")
print(distances(g, v = "U", to = "Y"))

library(causaleffect)
library(igraph)

# Step 1: Create the basic graph structure
g <- graph.formula(A -+ W -+ S, A -+ Y, U -+ W, U -+ Y)

# Step 2: Verify it's a DAG
cat("Is DAG:", is_dag(g), "\n")

# Step 3: Mark latent variable edges
g <- set_edge_attr(graph = g, 
                   name = "description", 
                   index = E(g)[.from("U")], 
                   value = "U")

# Step 4: Check edge attributes
cat("Edges from U:\n")
print(E(g)[.from("U")])
cat("Edge descriptions:\n")
print(edge_attr(g, "description"))

# Step 5: Try causal effect
result <- tryCatch({
  causal.effect("Y", "A", "S", G = g, steps = TRUE)
}, error = function(e) {
  cat("Error:", e$message, "\n")
  return(NULL)
})

if (!is.null(result)) {
  print("Success! Result:")
  print(result)
} else {
  cat("Failed - trying alternative approach\n")
}