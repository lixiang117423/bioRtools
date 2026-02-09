# Test script for plot_synteny function
# This script tests the plot_synteny function with example data

# Load the package
devtools::load_all()

# Test 1: Basic usage with default parameters
cat("=== Test 1: Basic usage ===\n")
result1 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link
)
print(result1$plot.synteny)
print(result1$data.summary)

# Test 2: Custom species labels
cat("\n=== Test 2: Custom species labels ===\n")
result2 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("Species A", "Species B", "Species C", "Species D")
)
print(result2$plot.synteny)

# Test 3: Custom colors
cat("\n=== Test 3: Custom colors ===\n")
result3 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D"),
  fill_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")
)
print(result3$plot.synteny)

# Test 4: With title and subtitle
cat("\n=== Test 4: With title and subtitle ===\n")
result4 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D"),
  title = "Gene Synteny Visualization",
  subtitle = "Example from bioRtools package"
)
print(result4$plot.synteny)

# Test 5: Hide gene labels
cat("\n=== Test 5: Hide gene labels ===\n")
result5 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D"),
  show_gene_labels = FALSE
)
print(result5$plot.synteny)

# Test 6: Adjust link transparency
cat("\n=== Test 6: Adjust link transparency ===\n")
result6 <- plot_synteny(
  gene_data = df.synteny.gene,
  syntenic_data = df.synteny.link,
  species_labels = c("A", "B", "C", "D"),
  link_alpha = 0.8
)
print(result6$plot.synteny)

cat("\n=== All tests completed successfully! ===\n")
