#' @export
#' @import ggplot2
#' @importFrom clusterSim data.Normalization
#' @import netbenchmark
#' @import gridExtra
#' @import data.table
#' @import glmnet
#' @import vegan
#' @import qgraph
#' @importFrom PANR assoScore
#'
#'@export
#'@name GRN
#'@alias GRN
#'@alias number_of_features
#'@alias support
#'@alias confidence
#'@title Gene regulatory network of cell morphological phenotypes
#'@description We use the mining association model to detect the associations between the landmark genes, and inference of gene regulatory network of cell morphological phenotypes. 
#'@author {Isar Nassiri, Matthew McCall}
#'@param number_of_features 
#'The number of top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis (e.g. 20).
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@param support
#'We use an association mining model to detect the associations between the genes and to infer a gene regulatory network from phenotypic experimental data. In order to select significant interactions between any two genes from
#'the complete digraph of interactions between genes, we use minimum thresholds on support and confidence measures. Support measures the frequency of an indicated gene in the cell morphological gene sets.
#'@param confidence 
#'We use an association mining model to detect the associations between the genes and to infer a gene regulatory network from phenotypic experimental data. In order to select significant interactions between any two genes from
#'the complete digraph of interactions between genes, we use minimum thresholds on support and confidence measures. Confidence shows how often a given association
#'rule between two genes has been found in the dataset.
#'@arguments {
#'@item{A Number}(number_of_features)
#'@item{A Number}(support)
#'@item{A Number}(confidence)
#'}
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'GRN(10, 0.1, 0.6)
#'@export

GRN <- function(number_of_features, support, confidence)
{
  L1000_TP_profiles <- Transcriptomic_Profile
  L1000_MP_profiles <- Cell_Morphology_Profile
  x_new <- Query_Transcriptomic_Profile
  
  L1000_TP_profiles <- scale(L1000_TP_profiles)
  L1000_MP_profiles <- scale(L1000_MP_profiles)
  
  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)
  
  performance <- data.frame()
  TP <- TN <- FP <- FN <- 0
  
  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a <- query_binary[i]  #my standard
      b <- repositoyr_binary[j,i]
      
      if(a+b == 2) { TP <- TP + 1 }
      if(a+b == 0) { TN <- TN + 1 }
      if(a == 0 && b == 1) { FP <- FP + 1 }
      if(a == 1 && b == 0) { FN <- FN + 1 }
    }
    performance[j,1] <- TP
    performance[j,2] <- TN
    performance[j,3] <- FP
    performance[j,4] <- FN
    performance[j,5] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    TP <- TN <- FP <- FN <- 0
  }
  
  colnames(performance) <- c("TP", "TN", "FP", "FN", "MCC")
  rownames(performance) <- rownames(repositoyr_binary)
  
  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)
  
  CMP_subset <<- L1000_MP_profiles[selected_drugs,]
  TP_subset <<- L1000_TP_profiles[selected_drugs,]
  dim(TP_subset)
  
  CMP_subset <- CMP_subset[-which(rownames(CMP_subset) %in% "BRD-K37798499"),]
  TP_subset <- TP_subset[-which(rownames(TP_subset) %in% "BRD-K37798499"),]
  
  MCR <- list();
  
  for(i in 1:number_of_features)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]
    
    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response
    
    set.seed(1)
    fit.lasso <- glmnet(x,y, standardize=TRUE)
    cv.lasso<-cv.glmnet(x,y)
    lam.best <- cv.lasso$lambda.min
    
    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=cv.lasso$lambda.min)
    Lasso_coefficient <- as.matrix(Lasso_coefficient)
    
    inds <- which(Lasso_coefficient[,1]!=0)
    variables<-row.names(Lasso_coefficient)[inds]
    variables<-variables[!(variables %in% '(Intercept)')];
    
    length(which(0 != Lasso_coefficient[,1]))
    results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
    colnames(results) <- "coef"
    results <- as.data.frame(results[-1,,FALSE])
    dim(results)
    selected_genes <- rownames(results)
    
    ste_name <- data.frame();
    
    LL <- NULL
    LL <- (length(selected_genes))
    
    if(0<length(selected_genes))
    {
      selected_GO <- selected_genes
      
      for(j2 in 1:LL)
      {
        ste_name[j2,1] <- selected_GO[j2]
      }
      
      if(!is.null(ste_name$V1))
      {
        MCR[[i]] <- ste_name[complete.cases(ste_name),1]
        names(MCR)[i] <- colnames(x_new2)
      } else {
        MCR[[i]] <- "NULL"
        names(MCR)[i] <- colnames(x_new2)
      }
    }
  }
  
  #-- Remove cm terms without any associated gene set
  null_cms <- which(sapply(MCR, is.null))
  if(0<length(null_cms)){ MCR <- MCR[-c(null_cms)] }
  #--
  
  df4 <- data.frame(c("gene","pathway"))
  
  for(i in 1:length(MCR))
  {
    df <- data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE)
    df2 <- rep(names(MCR[i]), dim(df)[2])
    df3 <-  rbind(df,df2)
    df4 <-  cbind(df3,df4)
  }
  
  length(names(MCR))
  
  df4 <- t(df4)
  colnames(df4) <- c("name","feature")
  df4 <- (df4[-dim(df4)[1],])
  dim(df4)
  
  #-- Association analysis
  
  all_genes_in_MCR <- as.data.frame(unique(df4[,1]))
  colnames(all_genes_in_MCR) <- "gene_name"
  dim(all_genes_in_MCR)
  
  #-- Support(X)
  
  in_gene_set <- 0
  
  for(j in 1:length(all_genes_in_MCR$gene_name))
  {
    
    for(i in 1:length(MCR))
    {
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE))
      class(df)
      if(0<length(intersect(as.character(all_genes_in_MCR$gene_name[j]), df))) { in_gene_set <- in_gene_set + 1  }
    }
    
    all_genes_in_MCR[j,2] <- in_gene_set/length(MCR)
    in_gene_set <- 0
  }
  
  all_genes_in_MCR_t <- all_genes_in_MCR[which(all_genes_in_MCR$V2 > support),]
  dim(all_genes_in_MCR_t)
  
  #-- Combination of genes --
  
  combination_all_genes_in_MCR <- expand.grid(all_genes_in_MCR_t$gene_name, all_genes_in_MCR_t$gene_name)
  class(combination_all_genes_in_MCR)
  colnames(combination_all_genes_in_MCR) <- c("gene1","gene2")
  
  #-- Remove self-loops --
  
  o <- 1
  
  selfloop <- data.frame()
  
  for(j in 1:dim(combination_all_genes_in_MCR)[1])
  {
    if(0 < length(which(as.character(combination_all_genes_in_MCR$gene1[j]) %in% as.character(combination_all_genes_in_MCR$gene2[j]))))
    { selfloop[o,1] <- j; o = o + 1 }
  }
  
  combination_all_genes_in_MCR <- combination_all_genes_in_MCR[-as.numeric(selfloop$V1),]
  dim(combination_all_genes_in_MCR)
  
  #--- Support(X,Y) + confidence ---
  
  in_gene_set <- 0
  
  for(j in 1:dim(combination_all_genes_in_MCR)[1])
  {
    
    for(i in 1:length(MCR))
    {
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE))
      
      if(2 == length(intersect(c(as.character(combination_all_genes_in_MCR$gene1[j]), as.character(combination_all_genes_in_MCR$gene2[j])),df)))
      { in_gene_set <- in_gene_set + 1  }
    }
    
    combination_all_genes_in_MCR[j,3] <- in_gene_set/length(MCR)
    
    combination_all_genes_in_MCR[j,4] <-
      (in_gene_set/length(MCR))/all_genes_in_MCR_t[which(all_genes_in_MCR_t$gene_name %in% as.character(combination_all_genes_in_MCR$gene1[j])),2]
    
    in_gene_set <- 0
  }
  
  colnames(combination_all_genes_in_MCR) <- c("gene1","gene2","supp(XvY)","Confidence")
  
  combination_all_genes_in_MCR_sig <- combination_all_genes_in_MCR[combination_all_genes_in_MCR$Confidence > confidence,]
  edgelist2 <- paste(combination_all_genes_in_MCR_sig[,1], combination_all_genes_in_MCR_sig[,2], sep = "~")
  length(edgelist2)
  
  # correlation network:
  clr <- clr.wrap(TP_subset)
  
  graph <- graph.adjacency(clr)
  edgelist <- get.edgelist(graph)
  edgelist1 <- paste(edgelist[,1], edgelist[,2], sep = "~")
  length(edgelist1)
  
  #-- comparison of rules and correlation
  
  selected_edges <- edgelist1[which(edgelist1 %in% edgelist2)]
  length(selected_edges)
  
  #Visulization of graph
  graph_e <- data.frame()
  
  for (i in 1:length(selected_edges))
  {
    e_ <- unlist(strsplit(selected_edges[i], "~"))
    graph_e[i,1] <- e_[1]
    graph_e[i,2] <- e_[2]
  }
  
  matrix_of_interactions <- as.matrix(graph_e)
  colnames(matrix_of_interactions) <- c("geneSymbol", "geneSymbol2")
  GeneRegulatoryNetwork <<- matrix_of_interactions[!duplicated(matrix_of_interactions),]
  dim(GeneRegulatoryNetwork)
  length(unique(c(graph_e[,1],graph_e[,2])))
  
  # -- Save data --
  return(GeneRegulatoryNetwork)
  
  print("You can find the results in the 'GeneRegulatoryNetwork' R object.")

}

#'@export
#'@name CellMorphologyEnrichmentAnalysis
#'@alias CellMorphologyEnrichmentAnalysis
#'@alias number_of_features
#'@title Cell morphology enrichment analysis
#'@description We use a stepwise variable selection approach, including combination of least absolute shrinkage and selection operator (LASSO) with cross-validation to tune parameter for cell morphology enrichment analysis. We consider all transcriptomic profiles in the reference repository as inputs of LASSO to select a subset of landmark genes (v) that best describe an indicated profile of cell morphology feature.
#'@author {Isar Nassiri, Matthew McCall}
#'@param number_of_features
#'The number of top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis (e.g. 20).
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@arguments item{A Number}(number_of_features)
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'CellMorphologyEnrichmentAnalysis(20)
#'@export

CellMorphologyEnrichmentAnalysis <- NULL
CellMorphologyEnrichmentAnalysis <- function(number_of_features)
{

  L1000_TP_profiles <- Transcriptomic_Profile
  L1000_MP_profiles <- Cell_Morphology_Profile
  x_new <- Query_Transcriptomic_Profile
  
  L1000_TP_profiles <- scale(L1000_TP_profiles)
  L1000_MP_profiles <- scale(L1000_MP_profiles)

  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)

  performance <- data.frame()
  TP <- TN <- FP <- FN <- 0

  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a <- query_binary[i]  #my standard
      b <- repositoyr_binary[j,i]

      if(a+b == 2) { TP <- TP + 1 }
      if(a+b == 0) { TN <- TN + 1 }
      if(a == 0 && b == 1) { FP <- FP + 1 }
      if(a == 1 && b == 0) { FN <- FN + 1 }
    }
    performance[j,1] <- TP
    performance[j,2] <- TN
    performance[j,3] <- FP
    performance[j,4] <- FN
    performance[j,5] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    TP <- TN <- FP <- FN <- 0
  }

  colnames(performance) <- c("TP", "TN", "FP", "FN", "MCC")
  rownames(performance) <- rownames(repositoyr_binary)

  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)

  CMP_subset <<- L1000_MP_profiles[selected_drugs,]
  TP_subset <<- L1000_TP_profiles[selected_drugs,]

  MCR <- list();

  for(i in 1:number_of_features)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]

    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response

    set.seed(1)
    fit.lasso <- glmnet(x,y, standardize=TRUE)
    cv.lasso<-cv.glmnet(x,y)
    lam.best <- cv.lasso$lambda.min

    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=cv.lasso$lambda.min)
    Lasso_coefficient <- as.matrix(Lasso_coefficient)

    inds <- which(Lasso_coefficient[,1]!=0)
    variables<-row.names(Lasso_coefficient)[inds]
    variables<-variables[!(variables %in% '(Intercept)')];

    length(which(0 != Lasso_coefficient[,1]))
    results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
    colnames(results) <- "coef"
    results <- as.data.frame(results[-1,,FALSE])
    dim(results)
    selected_genes <- rownames(results)

    ste_name <- data.frame();

    LL <- NULL
    LL <- (length(selected_genes))

    if(0<length(selected_genes))
    {
      selected_GO <- selected_genes

      for(j2 in 1:LL)
      {
        ste_name[j2,1] <- selected_GO[j2]
      }

      if(!is.null(ste_name$V1))
      {
        MCR[[i]] <- ste_name[complete.cases(ste_name),1]
        names(MCR)[i] <- colnames(x_new2)
      } else {
        MCR[[i]] <- "NULL"
        names(MCR)[i] <- colnames(x_new2)
      }
    }
  }

  null_cms <- which(sapply(MCR, is.null))
  if(0<length(null_cms)){ MCR <- MCR[-c(null_cms)] }
  #--

  df4 <- data.frame(c("gene","pathway"))

  for(i in 1:length(MCR))
  {
    df <- data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE)
    df2 <- rep(names(MCR[i]), dim(df)[2])
    df3 <-  rbind(df,df2)
    df4 <-  cbind(df3,df4)
  }

  length(names(MCR))

  df4 <- t(df4)
  colnames(df4) <- c("name","feature")
  datCM2 <- (df4[-dim(df4)[1],])

  agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))
  Names <- as.character(agregatation$feature)

  setDT(agregatation)[, id := .GRP, by = name]
  agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]
  agregatation$id <- sprintf("Cluster_%d", agregatation$id)
  CellMorphologyEnrichment <<- as.data.frame(agregatation)
  length(CellMorphologyEnrichment$feature)
  
  #-- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(CellMorphologyEnrichment)

  print("You can find the results in the 'CellMorphologyEnrichment' R objects.")

 }

#'@export
#'@name RankCellMorphologicalFeatures
#'@alias RankCellMorphologicalFeatures
#'@alias TOP
#'@title {Ranking of cell morphological phenotypes based on the Strength Centrality Score (SCS).}
#'@description We consider the connectivity and centrality of cell morphological features in Cosine similarity network of image-based cell morphological profile to rank them based on the Strength Centrality Score (SCS).
#'@author{Isar Nassiri, Matthew McCall}
#'@param TOP
#'The number of top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis (e.g. 20).
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@arguments item{A Number}(TOP)
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'RankCellMorphologicalFeatures(20)
#'@export

RankCellMorphologicalFeatures <- function(TOP)
{
  L1000_TP_profiles <- Transcriptomic_Profile
  L1000_MP_profiles <- Cell_Morphology_Profile
  x_new <- Query_Transcriptomic_Profile
  
  L1000_TP_profiles <- scale(L1000_TP_profiles)
  L1000_MP_profiles <- scale(L1000_MP_profiles)

  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)

  performance <- data.frame()
  TP <- TN <- FP <- FN <- 0

  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a <- query_binary[i]  #my standard
      b <- repositoyr_binary[j,i]

      if(a+b == 2) { TP <- TP + 1 }
      if(a+b == 0) { TN <- TN + 1 }
      if(a == 0 && b == 1) { FP <- FP + 1 }
      if(a == 1 && b == 0) { FN <- FN + 1 }
    }
    performance[j,1] <- TP
    performance[j,2] <- TN
    performance[j,3] <- FP
    performance[j,4] <- FN
    performance[j,5] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    TP <- TN <- FP <- FN <- 0
  }

  colnames(performance) <- c("TP", "TN", "FP", "FN", "MCC")
  rownames(performance) <- rownames(repositoyr_binary)

  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)

  CMP_subset <<- L1000_MP_profiles[selected_drugs,]
  TP_subset <<- L1000_TP_profiles[selected_drugs,]

  MCR <- list();

  for(i in 1:TOP)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]

    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response

    set.seed(1)
    fit.lasso <- glmnet(x,y, standardize=TRUE)
    cv.lasso<-cv.glmnet(x,y)
    lam.best <- cv.lasso$lambda.min

    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=cv.lasso$lambda.min)
    Lasso_coefficient <- as.matrix(Lasso_coefficient)

    inds <- which(Lasso_coefficient[,1]!=0)
    variables<-row.names(Lasso_coefficient)[inds]
    variables<-variables[!(variables %in% '(Intercept)')];

    length(which(0 != Lasso_coefficient[,1]))
    results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
    colnames(results) <- "coef"
    results <- as.data.frame(results[-1,,FALSE])
    dim(results)
    selected_genes <- rownames(results)

    ste_name <- data.frame();

    LL <- NULL
    LL <- (length(selected_genes))

    if(0<length(selected_genes))
    {
      selected_GO <- selected_genes

      for(j2 in 1:LL)
      {
        ste_name[j2,1] <- selected_GO[j2]
      }

      if(!is.null(ste_name$V1))
      {
        MCR[[i]] <- ste_name[complete.cases(ste_name),1]
        names(MCR)[i] <- colnames(x_new2)
      } else {
        MCR[[i]] <- "NULL"
        names(MCR)[i] <- colnames(x_new2)
      }
    }
  }

  null_cms <- which(sapply(MCR, is.null))
  if(0<length(null_cms)){ MCR <- MCR[-c(null_cms)] }
  #--

  df4 <- data.frame(c("gene","pathway"))

  for(i in 1:length(MCR))
  {
    df <- data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE)
    df2 <- rep(names(MCR[i]), dim(df)[2])
    df3 <-  rbind(df,df2)
    df4 <-  cbind(df3,df4)
  }

  length(names(MCR))

  df4 <- t(df4)
  colnames(df4) <- c("name","feature")
  datCM2 <- (df4[-dim(df4)[1],])

  agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))
  Names <- as.character(agregatation$feature)

  setDT(agregatation)[, id := .GRP, by = name]
  agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]   #sort based on the indices
  agregatation$id <- sprintf("Cluster_%d", agregatation$id) 			#add term of cluster at the beginning of each index
  agregatation <- as.data.frame(agregatation)
  length(agregatation$feature)

  #---------------------------
  CMP_subset2 <- as.matrix(CMP_subset[,which(colnames(CMP_subset) %in% agregatation$feature)])

  cor_bfi <- assoScore(CMP_subset2, "cosine", upperTri=FALSE, transform=FALSE)
  as <- centralityTable(cor_bfi)
  as_strength <- as[which(as$measure == "Strength"), ]

  colnames(as_strength) <- c("graph", "type", "feature", "variable", "Strength_centrality")
  as_strength <- as_strength[ , c(3,5)]
  agregatation <- merge(agregatation, as_strength, by = "feature")

  Long <- as_strength[order(as_strength[, 2], decreasing=TRUE)[1:TOP],]
  Long$type <- rep(1, dim(Long)[1])

  Long$feature <- factor(Long$feature, levels = Long$feature[order(Long$Strength_centrality)]) #I order y axis based on the x axis
  RankCellMorphological <<- Long[!is.na(Long$feature),]
  g <- ggplot(RankCellMorphological, aes(x = Strength_centrality, y = feature, group = type))
  g + xlab("") + ylab("") + geom_point() +  geom_path() + labs(list(title = "Top enriched cell morphological phenotypes", x = "Single-cell morphological phenotype score (Strength)", y = "Single-cell morphological phenotype")) + theme_grey(base_size = 10)

  #-- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(RankCellMorphological)
  
  print("You can find the results in the 'RankCellMorphological' R object.")

 }
  
#'@export   
#'@name crossTabulation
#'@alias crossTabulation
#'@alias TOP
#'@title {Cross-tabulation of landmark genes, and single-cell morphological features}
#'@description We present the results of cell morphology enrichment analysis as cross-tabulation of landmark genes, and single-cell morphological features including the direction of effects (up or down regulation).
#'@author {Isar Nassiri, Matthew McCall}
#'@param TOP
#'We use the top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis.
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@arguments {
#'@item{A Number}{
#'We use the top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis.
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'}(TOP)
#'}
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'crossTabulation(20) 
#'@export 
 
crossTabulation <- function(TOP)
{
L1000_TP_profiles <- Transcriptomic_Profile
L1000_MP_profiles <- Cell_Morphology_Profile
x_new <- Query_Transcriptomic_Profile

L1000_TP_profiles <- scale(L1000_TP_profiles)
L1000_MP_profiles <- scale(L1000_MP_profiles)

query_binary <- ifelse(x_new > 0, 1, 0)
repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)

performance <- data.frame()
TP <- TN <- FP <- FN <- 0

for(j in 1:dim(repositoyr_binary)[1])
{
  for(i in 1:dim(repositoyr_binary)[2])
  {
    a <- query_binary[i]  #my standard
    b <- repositoyr_binary[j,i]

    if(a+b == 2) { TP <- TP + 1 }
    if(a+b == 0) { TN <- TN + 1 }
    if(a == 0 && b == 1) { FP <- FP + 1 }
    if(a == 1 && b == 0) { FN <- FN + 1 }
  }
  performance[j,1] <- TP
  performance[j,2] <- TN
  performance[j,3] <- FP
  performance[j,4] <- FN
  performance[j,5] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  TP <- TN <- FP <- FN <- 0
}

colnames(performance) <- c("TP", "TN", "FP", "FN", "MCC")
rownames(performance) <- rownames(repositoyr_binary)

selected_drugs <- which(performance[,5] > 0.1)
length(selected_drugs)

CMP_subset <<- L1000_MP_profiles[selected_drugs,]
TP_subset <<- L1000_TP_profiles[selected_drugs,]

MCR <- list();

for(i in 1:TOP)
{
  x_new2 <- as.data.frame(CMP_subset[,i])
  colnames(x_new2) <- colnames(CMP_subset)[i]

  x <- data.matrix(TP_subset)                #predictors
  y <- as.numeric(unlist(x_new2))            #response

  set.seed(1)
  fit.lasso <- glmnet(x,y, standardize=TRUE)
  cv.lasso<-cv.glmnet(x,y)
  lam.best <- cv.lasso$lambda.min

  #-- extract significant coefficients --
  Lasso_coefficient <- coef(fit.lasso, s=cv.lasso$lambda.min)
  Lasso_coefficient <- as.matrix(Lasso_coefficient)

  inds <- which(Lasso_coefficient[,1]!=0)
  variables<-row.names(Lasso_coefficient)[inds]
  variables<-variables[!(variables %in% '(Intercept)')];

  length(which(0 != Lasso_coefficient[,1]))
  results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
  colnames(results) <- "coef"
  results <- as.data.frame(results[-1,,FALSE])
  dim(results)
  selected_genes <- rownames(results)

  ste_name <- data.frame();

  LL <- NULL
  LL <- (length(selected_genes))

  if(0<length(selected_genes))
  {
    selected_GO <- selected_genes

    for(j2 in 1:LL)
    {
      ste_name[j2,1] <- selected_GO[j2]
    }

    if(!is.null(ste_name$V1))
    {
      MCR[[i]] <- ste_name[complete.cases(ste_name),1]
      names(MCR)[i] <- colnames(x_new2)
    } else {
      MCR[[i]] <- "NULL"
      names(MCR)[i] <- colnames(x_new2)
    }
  }
}

null_cms <- which(sapply(MCR, is.null))
if(0 < length(null_cms)){ MCR <- MCR[-c(null_cms)] }
#--

df4 <- data.frame(c("gene","pathway"))

for(i in 1:length(MCR))
{
  df <- data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE)
  df2 <- rep(names(MCR[i]), dim(df)[2])
  df3 <-  rbind(df,df2)
  df4 <-  cbind(df3,df4)
}

length(names(MCR))

df4 <- t(df4)
colnames(df4) <- c("name","feature")
datCM2 <- (df4[-dim(df4)[1],])

agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))
Names <- as.character(agregatation$feature)

setDT(agregatation)[, id := .GRP, by = name]   						#Package 'data.table' - CRAN, it indexes the gene sets (add new column of indices)
agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]   #Sort based on the indices
agregatation$id <- sprintf("Cluster_%d", agregatation$id)             #Add term of cluster at the beginning of each index
agregatation <- as.data.frame(agregatation)
length(agregatation$feature)

#---------------------------
CMP_subset2 <- as.matrix(CMP_subset[,which(colnames(CMP_subset) %in% agregatation$feature)])

cor_bfi <- assoScore(CMP_subset2, "cosine", upperTri=FALSE, transform=FALSE)
as <- centralityTable(cor_bfi)
as_strength <- as[which(as$measure == "Strength"), ]

colnames(as_strength) <- c("graph", "type", "feature", "variable", "Strength_centrality")
as_strength <- as_strength[ , c(3,5)]
agregatation <- merge(agregatation, as_strength, by = "feature")

Long <- as_strength[order(as_strength[, 2], decreasing=TRUE)[1:TOP],]
Long$type <- rep(1, dim(Long)[1])

#-- visulization of crosstab table ---
datCM2 <- as.data.frame(datCM2)
df5 <- datCM2[which(as.character(datCM2$feature) %in% as.character(Long$feature)),]
df5 <- as.data.frame(df5)
length(unique(as.character(df5[,2])))

query <- as.data.frame(Query_Transcriptomic_Profile)

for(i in 1:dim(df5)[1])
{
  df5[i,3] <- query[which(rownames(query) %in% df5[i,1]),]
}

colnames(df5) <- c("Gene_name","Cell_morphological_phenotype","weight")

weights <- as.numeric(df5$weight)

for(i in 1:length(weights))
{
  if(weights[i] > 0) {weights[i] <- 1}
  if(weights[i] < 0) {weights[i] <- -1}
}

df5$weight <- weights

ggplot(df5, aes(x = Gene_name, y = Cell_morphological_phenotype, fill = weight)) +
  geom_raster() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 270, hjust = 0),
        axis.text=element_text(size=8), axis.title=element_text(size=14, face="bold"))

df6 <- data.frame()

for (num_row in 1:dim(df5)[1])
{
  df6[num_row,1] <- as.character(levels(df5[num_row,1])[df5[num_row,1]])
  df6[num_row,2] <- as.character(levels(df5[num_row,2])[df5[num_row,2]])
  df6[num_row,3] <- df5[num_row,3]
}

M1 <- as.matrix(df6)
G1 <- graph.data.frame(M1)
MA <- get.adjacency(G1,sparse=FALSE)

MA[M1[,1:2]] <- as.numeric(M1[,3])  #weights

crossTable <<- MA[1:length(unique(df6[,1])),-(1:length(unique(df6[,1])))]

#-- Save data --
return(CMP_subset)
return(TP_subset)
return(crossTable)
  
print("You can find the results in the 'crossTable' R object")

}

#'@name Mapping
#'@alias Mapping
#'@title Mapping query transcriptomic profile against the reference repository
#'@description {We map the query profile (Q) of fold change of gene expression versus reference repository of transcriptomic data (T) to detect similarities among these profiles (s).
#'First, we convert the query and backend repository of transcription profiles to the Boolean expression and replace up-regulated values with one and down regulated
#'values with zero. We provide a confusion matrix for the intersection of the query with each of the transcriptomic profiles in the reference repository.
#'Next, we use the confusion matrices to score the similarity between the query and reference transcription profiles based on the Matthew correlation value.
#'}
#'@author {Isar Nassiri, Matthew McCall}
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'Mapping()
#'@export

Mapping <- NULL
Mapping <- function()
{

  L1000_TP_profiles <- Transcriptomic_Profile
  L1000_MP_profiles <- Cell_Morphology_Profile
  x_new <- Query_Transcriptomic_Profile
  
  L1000_TP_profiles <- scale(L1000_TP_profiles)
  L1000_MP_profiles <- scale(L1000_MP_profiles)

  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)

  performance <- data.frame()
  TP <- TN <- FP <- FN <- 0

  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a <- query_binary[i]  #my standard
      b <- repositoyr_binary[j,i]

      if(a+b == 2) { TP <- TP + 1 }
      if(a+b == 0) { TN <- TN + 1 }
      if(a == 0 && b == 1) { FP <- FP + 1 }
      if(a == 1 && b == 0) { FN <- FN + 1 }
    }
    performance[j,1] <- TP
    performance[j,2] <- TN
    performance[j,3] <- FP
    performance[j,4] <- FN
    performance[j,5] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    TP <- TN <- FP <- FN <- 0
  }

  colnames(performance) <- c("TP", "TN", "FP", "FN", "MCC")
  rownames(performance) <- rownames(repositoyr_binary)

  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)

  CMP_subset <<- L1000_MP_profiles[selected_drugs,]
  TP_subset <<- L1000_TP_profiles[selected_drugs,]

  #-- Save data --
  return(CMP_subset)
  return(TP_subset)

  print("You can find the results in the CMP_subset and TP_subset R objects.")

}