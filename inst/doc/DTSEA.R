## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DTSEA)

## ----echo = FALSE-------------------------------------------------------------
# Load the data
data("example_disease_list", package = "DTSEA")
data("example_drug_target_list", package = "DTSEA")
data("example_ppi", package = "DTSEA")

# Perform a simple DTSEA analysis using default optional parameters. 
result <- DTSEA(network = example_ppi,
                disease = example_disease_list,
                drugs = example_drug_target_list
)

## -----------------------------------------------------------------------------
library(dplyr)
select(result, -leadingEdge) %>%
  arrange(desc(NES)) %>%
  filter(NES > 0 & pval < .05)

## -----------------------------------------------------------------------------
# Calculate p0
p0 <- calculate_p0(nodes = example_ppi, disease = example_disease_list)

# Then perform random walk
random.walk(network = example_ppi, p0 = p0)

## ----echo = FALSE-------------------------------------------------------------
# Imagine there are three prediction results with ten samples
x <- runif(10, min = 0, max = 5)
y <- runif(10, min = 0, max = 5)
z <- sqrt(x + y) + runif(10, min = -1, max = 1)
data <- data.frame(x, y, z)

## -----------------------------------------------------------------------------
# Just report the results
kendall.w(data)$report

# Or just report the alpha
cronbach.alpha(data)

