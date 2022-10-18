#' @title Main function of drug target set enrichment analysis (DTSEA)
#' @description The DTSEA function determines whether a drug is potent for a
#' specific disease by the proximity between its targets and the disease-related
#' genes.
#' @param network The human protein-protein interactome network. It should be or
#' be preconverted before being inputted in DTSEA.
#' @param disease The disease-related nodes.
#' @param drugs The drug-target long format dataframe. It includes at least
#' columns with the drug_id and drug_target.
#' @param rwr.pt The random walk p0 vector. Set it to 0 if you wish DTSEA
#' automatically compute it, or you can provide your predetermined p0 vector.
#' @param eps The boundary of calculating the p value.
#' @param nPermSimple Number of permutations in the simple fgsea implementation
#' for preliminary estimation of P-values.
#' @param gseaParam GSEA parameter value, all gene-level statistics are raised
#' to the power of 'gseaParam' before calculating of GSEA enrichment scores.
#' @param verbose Show the messages
#'
#' @importFrom fgsea fgseaMultilevel
#' @importFrom igraph graph_from_data_frame graph_from_edgelist simplify induced_subgraph
#' @importFrom dplyr rename %>%
#' @importFrom tibble as_tibble
#'
#' @return The resulting dataframe consists of `drug_id`, `pval`, `padj`,
#' `log2err`, `ES`, `NES`, `size`, and `leadingEdge`.
#' @export
#'
#' @examples
#' library(dplyr)
#' library(DTSEA)
#'
#' # Load the data
#' data("example_disease_list", package = "DTSEA")
#' data("example_drug_target_list", package = "DTSEA")
#' data("example_ppi", package = "DTSEA")
#'
#' # Run the DTSEA and sort the result dataframe by normalized enrichment scores
#' # (NES)
#' result <- DTSEA(
#'   network = example_ppi,
#'   disease = example_disease_list,
#'   drugs = example_drug_target_list
#' ) %>%
#' arrange(desc(NES))
#'
#'
#' # We can extract the significantly NES > 0 drug items.
#' result %>%
#'   filter(NES > 0 & pval < .05)
#' # Or we can draw the enrichment plot of the first predicted drug.
#' fgsea::plotEnrichment(
#'   pathway = example_drug_target_list %>%
#'     filter(drug_id == slice(result, 1)$drug_id) %>%
#'     pull(gene_target),
#'   stats = random.walk(network = example_ppi,
#'                       p0 = calculate_p0(nodes = example_ppi,
#'                                         disease = example_disease_list)
#'                       )
#' )
#'
DTSEA <- function(network, disease, drugs, rwr.pt = 0, eps = 1e-50,
                  nPermSimple = 5000, gseaParam = 1, verbose = TRUE) {

  # First, we need to simplify the whole graph
  if (is.data.frame(network)) {
    graph <- graph_from_data_frame(network, directed = FALSE)
  } else if (is.matrix(network)) {
    graph <- graph_from_edgelist(network, directed = FALSE)
  } else if (!is.igraph(network)) {
    stop("Not an igraph object")
  } else {
    graph <- network
  }
  graph <- simplify(graph, remove.loops = TRUE, remove.multiple = TRUE) %>%
    # Just want the main component
    induced_subgraph(V(.)[components(.)$membership == which.max(components(.)$csize)])

  drug.genes <- split(drugs$gene_target, drugs$drug_id)

  # Then we compute the p0 vector
  if (!(sum(rwr.pt))) {
    if (verbose) {
      message("Random walking...\n")
    }
    rwr.pt <- calculate_p0(nodes = graph, disease = disease) %>%
      random.walk(graph, ., verbose = verbose)
  }
  if (verbose) {
    message("Doing GSEA enrichment...\n")
  }
  # Do GSEA enrichment

  # The enrichment analysis aims to determine whether member S (drug genes)
  # are randomly distributed through L (RWR gene set).
  # Usually, S is biological pathways, while L is gene expression data.
  result <- fgseaMultilevel(
    pathways = drug.genes, stats = rwr.pt, eps = eps,
    nPermSimple = nPermSimple, gseaParam = gseaParam
  ) %>%
    as_tibble() %>%
    rename(drug_id = 1)

  return(result)
}
