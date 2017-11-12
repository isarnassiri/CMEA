#' @export
#' @importFrom netbenchmark clr.wrap
#' @import data.table
#' @importFrom glmnet glmnet
#' @importFrom reshape2 melt
#' @import igraph
#'
#'@export
#'@name geneInteractionNetwrok
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
#'@return 1. A matrix of gene-gene interaction network, including three columns: source gene, sink gene, and type of interaction (+ activation, - inhibition); 2. A data frame including profiles of cell morphological features of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column cell morphological features); 3. A data frame including gene expression profiles of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column gene symbol).
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'geneInteractionNetwrok(10, 0.1, 0.6)
#'@export

geneInteractionNetwrok <- NULL
geneInteractionNetwrok <- function(number_of_features, support, confidence)
{
  L1000_TP_profiles <- scale(Transcriptomic_Profile)
  L1000_MP_profiles <- scale(Cell_Morphology_Profile)

  x_new <- Query_Transcriptomic_Profile

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
   
  CMP_subset <- L1000_MP_profiles[selected_drugs,]
  TP_subset <- L1000_TP_profiles[selected_drugs,]
  
  #-- Data Normalization [unitization with zero minimum]
  for(l in 1:dim(CMP_subset)[2])
	{
	  min <- min(CMP_subset[,l])
	  max <- max(CMP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(CMP_subset)[1])
	  {  CMP_subset[t,l] <- ((CMP_subset[t,l]-min)/range) }
	}

  for(l in 1:dim(TP_subset)[2])
	{
	  min <- min(TP_subset[,l])
	  max <- max(TP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(TP_subset)[1])
	  {  TP_subset[t,l] <- ((TP_subset[t,l]-min)/range) }
	} 

  MCR <- list();

  for(i in 1:number_of_features)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]

    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response

    set.seed(1)
	fit.lasso <- glmnet(x, y, standardize=TRUE)

	#--- select Lambda ------
	train=sample(seq(dim(x)[1]),(dim(x)[1]/2),replace=FALSE)
	lasso.tr=glmnet(x[train,],y[train])
	pred=predict(lasso.tr,x[-train,])
	rmse= sqrt(apply((y[-train]-pred)^2,2,mean))
	#plot(log(lasso.tr$lambda),rmse,type="b",xlab="Log(lambda)")
	lam.best=lasso.tr$lambda[order(rmse)[1]]
	lam.best
	  
    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=lam.best)
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

  #-- Association analysis

  all_genes_in_MCR <- as.data.frame(unique(df4[,1]))
  colnames(all_genes_in_MCR) <- "gene_name"
  dim(all_genes_in_MCR)

  #-- n_x
  
  in_gene_set <- 0 
  
  for(j in 1:length(all_genes_in_MCR$gene_name))
  {
    
    for(i in 1:length(MCR))
    {
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE))
      class(df)
      if(0<length(intersect(as.character(all_genes_in_MCR$gene_name[j]), df))) { in_gene_set <- in_gene_set + 1  }
    }
    
    all_genes_in_MCR[j,2] <- in_gene_set
    in_gene_set <- 0 
  }

  dim(all_genes_in_MCR)
  
  all_genes_in_MCR <- all_genes_in_MCR[all_genes_in_MCR$V2/length(MCR) >= 0.05,]
  
  dim(all_genes_in_MCR)   #n_x
  
  #-- Combination of genes --
  
  combination_all_genes_in_MCR <- expand.grid(all_genes_in_MCR$gene_name, all_genes_in_MCR$gene_name)
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
  
  in_gene_set <- 0 
  
  for(j in 1:dim(combination_all_genes_in_MCR)[1])
  {
    
    for(i in 1:length(MCR))
    {
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=TRUE), stringsAsFactors=FALSE))
       
      if(2 == length(intersect(c(as.character(combination_all_genes_in_MCR$gene1[j]), as.character(combination_all_genes_in_MCR$gene2[j])),df)))
        { in_gene_set <- in_gene_set + 1  }
    }
    
    combination_all_genes_in_MCR[j,3] <- in_gene_set/length(MCR)  #n_ab/n
    
    combination_all_genes_in_MCR[j,4] <-
    in_gene_set/all_genes_in_MCR[which(all_genes_in_MCR$gene_name %in% as.character(combination_all_genes_in_MCR$gene1[j])),2]
    #n_ab/n_a
    combination_all_genes_in_MCR[j,5] <- combination_all_genes_in_MCR[j,3]/(all_genes_in_MCR$V2[which(all_genes_in_MCR$gene_name %in% as.character(combination_all_genes_in_MCR$gene1[j]))]/length(MCR)* 
    all_genes_in_MCR$V2[which(all_genes_in_MCR$gene_name %in% as.character(combination_all_genes_in_MCR$gene2[j]))]/length(MCR)) 
    
    in_gene_set <- 0 
  }
  
  colnames(combination_all_genes_in_MCR) <- c("gene1","gene2","supp(XvY)","Confidence", "lift")
  dim(combination_all_genes_in_MCR)
  
  # combination_all_genes_in_MCR_sig <- combination_all_genes_in_MCR[combination_all_genes_in_MCR$Confidence > 0.6,]
  # combination_all_genes_in_MCR_sig <- combination_all_genes_in_MCR_sig[combination_all_genes_in_MCR_sig$`supp(XvY)` > 0.1]
  # 
  combination_all_genes_in_MCR_sig <- combination_all_genes_in_MCR[combination_all_genes_in_MCR$lift > 4,]
  
  edgelist2 <- paste(combination_all_genes_in_MCR_sig[,1], combination_all_genes_in_MCR_sig[,2], sep = "~")
  length(edgelist2)
  
  # correlation network:
  clr <- clr.wrap(TP_subset)

  edgelist <- melt(clr)
  edgelist <- edgelist[which(edgelist$value != 0),]
  edgelist1 <- paste(edgelist[,1], edgelist[,2], sep = "~")

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
  matrix_of_interactions <- matrix_of_interactions[!duplicated(matrix_of_interactions),]
  
  # -- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(selected_edges)
  
  print("You can find the results in the 'GeneRegulatoryNetwork' R object.")

}

#'@export
#'@name cellMorphologyEnrichmentAnalysis
#'@title cell morphology enrichment analysis
#'@description We use a stepwise variable selection approach, including combination of least absolute shrinkage and selection operator (LASSO) with cross-validation to tune parameter for cell morphology enrichment analysis. We consider all transcriptomic profiles in the reference repository as inputs of LASSO to select a subset of landmark genes (v) that best describe an indicated profile of cell morphology feature.
#'@author {Isar Nassiri, Matthew McCall}
#'@param number_of_features
#'The number of top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis (e.g. 20).
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@return 1. A data frame including 812 cell morphological features (row) and associated land mark genes with each feature (column); 2. A data frame including profiles of cell morphological features of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column cell morphological features); 3. A data frame including gene expression profiles of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column gene symbol).
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'cellMorphologyEnrichmentAnalysis(20)
#'@export

cellMorphologyEnrichmentAnalysis <- NULL
cellMorphologyEnrichmentAnalysis <- function(number_of_features)
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
    
  CMP_subset <- L1000_MP_profiles[selected_drugs,]
  TP_subset <- L1000_TP_profiles[selected_drugs,]
    
  #-- Data Normalization [unitization with zero minimum]
  for(l in 1:dim(CMP_subset)[2])
	{
	  min <- min(CMP_subset[,l])
	  max <- max(CMP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(CMP_subset)[1])
	  {  CMP_subset[t,l] <- ((CMP_subset[t,l]-min)/range) }
	}

  for(l in 1:dim(TP_subset)[2])
	{
	  min <- min(TP_subset[,l])
	  max <- max(TP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(TP_subset)[1])
	  {  TP_subset[t,l] <- ((TP_subset[t,l]-min)/range) }
	} 

  MCR <- list();

  for(i in 1:number_of_features)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]

    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response

    set.seed(1)
	fit.lasso <- glmnet(x, y, standardize=TRUE)
	  
	#--- select Lambda ------
	train=sample(seq(dim(x)[1]),(dim(x)[1]/2),replace=FALSE)
	lasso.tr=glmnet(x[train,],y[train])
	pred=predict(lasso.tr,x[-train,])
	rmse= sqrt(apply((y[-train]-pred)^2,2,mean))
	#plot(log(lasso.tr$lambda),rmse,type="b",xlab="Log(lambda)")
	lam.best=lasso.tr$lambda[order(rmse)[1]]
	lam.best
	  
    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=lam.best)
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
  CellMorphologyEnrichment <- as.data.frame(agregatation)
  length(CellMorphologyEnrichment$feature)
  
  #-- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(CellMorphologyEnrichment)

  print("You can find the results in the 'CellMorphologyEnrichment' R objects.")

 }
  
#'@name mappingQueryTranscriptomic
#'@title Mapping query transcriptomic profile against the reference repository
#'@description {We map the query profile (Q) of fold change of gene expression versus reference repository of transcriptomic data (T) to detect similarities among these profiles (s).
#'First, we convert the query and backend repository of transcription profiles to the Boolean expression and replace up-regulated values with one and down regulated
#'values with zero. We provide a confusion matrix for the intersection of the query with each of the transcriptomic profiles in the reference repository.
#'Next, we use the confusion matrices to score the similarity between the query and reference transcription profiles based on the Matthew correlation value.
#'}
#'@author {Isar Nassiri, Matthew McCall}
#'@param Transcriptomic_Profile
#'A data frame including expression level of 978 land mark genes in response to treatment with 162 drugs/small compound molecules (row drug/small molecule compound ID, and column land mark gene symbols)
#'@param Cell_Morphology_Profile
#'A data frame including profiles of 812 cell morphological features of 162 drugs/small compound molecules (row drug/small molecule compound ID, and column cell morphological features)
#'@param Query_Transcriptomic_Profile
#'A data frame including expression level of 978 land mark genes in response to treatment with an indicated drugs/small compound molecule (row drug/small molecule compound ID, and column land mark gene symbols)
#'@return 1. A data frame including profiles of cell morphological features of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column cell morphological features); 3. A data frame including gene expression profiles of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column gene symbol).
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'mappingQueryTranscriptomic(Query_Transcriptomic_Profile, Transcriptomic_Profile, Cell_Morphology_Profile)
#'@export

mappingQueryTranscriptomic <- NULL
mappingQueryTranscriptomic <- function(Query_Transcriptomic_Profile, Transcriptomic_Profile, Cell_Morphology_Profile)
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
 
  CMP_subset <- L1000_MP_profiles[selected_drugs,]
  TP_subset <- L1000_TP_profiles[selected_drugs,]
  
  #-- Data Normalization [unitization with zero minimum]
  for(l in 1:dim(CMP_subset)[2])
	{
	  min <- min(CMP_subset[,l])
	  max <- max(CMP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(CMP_subset)[1])
	  {  CMP_subset[t,l] <- ((CMP_subset[t,l]-min)/range) }
	}

  for(l in 1:dim(TP_subset)[2])
	{
	  min <- min(TP_subset[,l])
	  max <- max(TP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(TP_subset)[1])
	  {  TP_subset[t,l] <- ((TP_subset[t,l]-min)/range) }
	} 

  #-- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(selected_drugs)

  print("You can find the results in R object under title of CMP_subset and TP_subset.")
}

#'@export   
#'@name crossTabulation
#'@title {Cross-tabulation of landmark genes, and single-cell morphological features}
#'@description We present the results of cell morphology enrichment analysis as cross-tabulation of landmark genes, and single-cell morphological features including the direction of effects (up or down regulation).
#'@author {Isar Nassiri, Matthew McCall}
#'@param TOP
#'We use the top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis.
#'This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.
#'@return 1. A crosstab matrix including 812 cell morphological features (row), associated land mark genes with each feature (column), and gene expression level vaules as its elements; 2. A data frame including profiles of cell morphological features of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column cell morphological features); 3. A data frame including gene expression profiles of similar drugs/small compound molecules with query (row drug/small molecule compound ID, and column gene symbol).
#'@examples
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'data(Query_Transcriptomic_Profile)
#'crossTabulation(20) 
#'@export 

crossTabulation <- NULL
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

CMP_subset <- L1000_MP_profiles[selected_drugs,]
TP_subset <- L1000_TP_profiles[selected_drugs,]

#-- Data Normalization [unitization with zero minimum]
  for(l in 1:dim(CMP_subset)[2])
	{
	  min <- min(CMP_subset[,l])
	  max <- max(CMP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(CMP_subset)[1])
	  {  CMP_subset[t,l] <- ((CMP_subset[t,l]-min)/range) }
	}

  for(l in 1:dim(TP_subset)[2])
	{
	  min <- min(TP_subset[,l])
	  max <- max(TP_subset[,l])
	  range <- max-min
	  
	  for(t in 1:dim(TP_subset)[1])
	  {  TP_subset[t,l] <- ((TP_subset[t,l]-min)/range) }
	} 

MCR <- list();

for(i in 1:TOP)
{
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]

    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response

    set.seed(1)
	fit.lasso <- glmnet(x, y, standardize=TRUE)
	  
	#--- select Lambda ------
	train=sample(seq(dim(x)[1]),(dim(x)[1]/2),replace=FALSE)
	lasso.tr=glmnet(x[train,],y[train])
	pred=predict(lasso.tr,x[-train,])
	rmse= sqrt(apply((y[-train]-pred)^2,2,mean))
	#plot(log(lasso.tr$lambda),rmse,type="b",xlab="Log(lambda)")
	lam.best=lasso.tr$lambda[order(rmse)[1]]
	lam.best
	  
    #-- extract significant coefficients --
    Lasso_coefficient <- coef(fit.lasso, s=lam.best)
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

setDT(agregatation)[, id := .GRP, by = name]   						      #Package 'data.table' - CRAN, it indexes the gene sets (add new column of indices)
agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]   #Sort based on the indices
agregatation$id <- sprintf("Cluster_%d", agregatation$id)                 #Add term of cluster at the beginning of each index
agregatation <- as.data.frame(agregatation)
length(agregatation$feature)

#-- visulization of crosstab table ---
df5 <- as.data.frame(datCM2)
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

#ggplot(df5, aes(x = Gene_name, y = Cell_morphological_phenotype, fill = weight)) +
#  geom_raster() +
#  theme_bw() +
#  theme(axis.text.x = element_text(angle = 270, hjust = 0),
#        axis.text=element_text(size=8), axis.title=element_text(size=14, face="bold"))

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

crossTable <- MA[1:length(unique(df6[,1])),-(1:length(unique(df6[,1])))]

#-- Save data --
return(CMP_subset)
return(TP_subset)
return(crossTable)
  
print("You can find the results in the 'crossTable' R object")

}
