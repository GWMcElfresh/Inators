
#' @title GWsFavoriteKEGGPathways
#' A function to return a vector of KEGG pathways that are commonly studied in GW's research.
#' @param KEGG_IDs Logical, if TRUE returns KEGG pathway IDs, if FALSE returns pathway names. Default is FALSE.
#'
#' @return A vector of KEGG pathway IDs or names.
#' @export
GWsFavoriteKEGGPathways <-function(KEGG_IDs = FALSE){
  if (KEGG_IDs) {
  KEGG_pathways <-  c("hsa04010", # MAPK signaling pathway
                      "hsa04110", # Cell cycle
                      "hsa04115", # p53 signaling pathway
                      "hsa04310", # Wnt signaling pathway
                      "hsa04350", # TGF-beta signaling pathway
                      "hsa04630", # JAK-STAT signaling pathway
                      "hsa04660", # T cell receptor signaling pathway
                      "hsa04670", # Leukocyte transendothelial migration
                      "hsa04722", # Neurotrophin signaling pathway
                      "hsa04910",  # Insulin signaling pathway
                      "hsa04931", # Viral protein interaction with cytokine and cytokine receptor
                      "hsa04060", # Cytokine-cytokine receptor interaction
                      "hsa04062", # Chemokine signaling pathway
                      "hsa04620", # Toll-like receptor signaling pathway
                      "hsa04621", # NOD-like receptor signaling pathway
                      "hsa04622", # RIG-I-like receptor signaling pathway
                      "hsa04064", # NF-kappa B signaling pathway
                      "hsa04630", # JAK-STAT signaling pathway
                      "hsa04210", # Apoptosis
                      "hsa04658", # Th1 and Th2 cell differentiation
                      "hsa04659", # Th17 cell differentiation
                      "hsa05205"  # PD-L1 expression and PD-1 checkpoint pathway in cancer
  )} else {
    KEGG_pathways <-  c("MAPK signaling pathway",
                        "Cell cycle",
                        "p53 signaling pathway",
                        "Wnt signaling pathway",
                        "TGF-beta signaling pathway",
                        "JAK-STAT signaling pathway",
                        "T cell receptor signaling pathway",
                        "Leukocyte transendothelial migration",
                        "Neurotrophin signaling pathway",
                        "Insulin signaling pathway",
                        "Viral protein interaction with cytokine and cytokine receptor",
                        "Cytokine-cytokine receptor interaction",
                        "Chemokine signaling pathway",
                        "Toll-like receptor signaling pathway",
                        "NOD-like receptor signaling pathway",
                        "RIG-I-like receptor signaling pathway",
                        "NF-kappa B signaling pathway",
                        "JAK-STAT signaling pathway",
                        "Apoptosis",
                        "Th1 and Th2 cell differentiation",
                        "Th17 cell differentiation",
                        "PD-L1 expression and PD-1 checkpoint pathway in cancer"
    )
  }
return(KEGG_pathways)
}
