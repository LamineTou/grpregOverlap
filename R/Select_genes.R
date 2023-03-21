#' Selction of the genes present in reactome pathways 
#'
#' @import data.table
#' @import dplyr
#' @import msigdbr
#' @param df Dataframe
#' @param ensembl2ggnc contain the ensembl gene names with their correspondance to hgnc gene symbol.
#'
#' @return genes included or not in the pathways as separated dataframe and groups as list for included genes.
select.genes <- function(ensembl2hgnc, df = NULL) {
  
  output <- list()
  # Extract pathways 
  reactome_gene_set <- Format_GS(Species = "Homo sapiens", cat = "C2", subcat = "CP:REACTOME")
  
  gens=names(df)
  genes1=ensembl2hgnc[ensembl2hgnc$V2 %in% gens,] 
  ####### Unlist pathways in the gene length
  uns=as.data.frame(with(reactome_gene_set, rep(gs_name, lengths(ensembl_gene))))

  ######## Unlist gene 
  unlisgen=as.data.frame(unlist(reactome_gene_set$ensembl_gene))
  
  rownames(uns)=rownames(unlisgen)
  ######### Merge gene with following pathways as dataframe
  pathways=merge(uns,unlisgen, by="row.names", all = TRUE)
  pathways=pathways[,-1]
  colnames(pathways)[which(names(pathways) == "with(reactome_gene_set, rep(gs_name, lengths(ensembl_gene)))")] <- "pathways"
  colnames(pathways)[which(names(pathways) == "unlist(reactome_gene_set$ensembl_gene)")] <- "genes"
  
  ####### Select only genes in df with hgnc symbol 
  gen=genes1$V1
  sub_path <- pathways[pathways$genes %in% gen,] 
  
  
  ### The list of sub set of genes present in the pathways
  gene=sub_path$genes
  sub_genes <- genes1[genes1$V1 %in% gene,]
  
  ### The list of sub set of genes absent in the pathways
  `%notin%` <- Negate(`%in%`)
  sub_genes_absent= genes1[genes1$V1 %notin% gene,]
  
  ######### Here is the data containing the gene present in the pathways. 
  subgen = sub_genes$V2
  gene_present=df[names(df) %in% subgen]
  
  ## Create the data of the set of absent gene from RNAseq data 
  subgen1 = sub_genes_absent$V2
  genes_abs <- df[names(df) %in% subgen1]
    
  output$include_genes <- gene_present
  output$absent_genes <- genes_abs
  
  # Regroup genes 
  gene_pathways_set <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
    dplyr::group_by(gs_name) %>%
    dplyr::select(gs_name, gene_symbol, human_ensembl_gene) %>%
    dplyr::mutate(pathways = list(gs_name)) %>%
    dplyr::select(gene_symbol, pathways, gs_name) %>%
    dplyr::distinct()
  
   
  gene_pathways_set = gene_pathways_set[,-2] 
  # subgen is the list of genes present in pathways
  my_group = gene_pathways_set[gene_pathways_set$gene_symbol %in% subgen,] 
  sub_pathway = my_group %>%
    group_by(gs_name) %>%
    summarise(pathways = list(gene_symbol))
  
  output$groups <- sub_pathway
  
  return(output)
}

# -------------------------------------------------------------------------------
## function: Format_GS select genes present in msigdb pathways
# -------------------------------------------------------------------------------

Format_GS <- function(Species = "Homo sapiens", 
                      cat = NULL, 
                      subcat = NULL) 
{
  require(msigdbr, quietly = T)
  require(dplyr, quietly = T)
  gene_set <- msigdbr(species = Species, category = cat, subcategory = subcat) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::select(gs_name, gene_symbol, human_ensembl_gene) %>%
    dplyr::mutate(ensembl_gene = strsplit(paste(human_ensembl_gene, collapse=","),",")) %>%
    dplyr::select(gs_name, ensembl_gene) %>%
    dplyr::distinct()
  return(gene_set)
}




