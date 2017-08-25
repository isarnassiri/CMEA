#' @export
#' @importFrom clusterSim data.Normalization
#' @importFrom netbenchmark clr.wrap
#' @importFrom data.table setDT
#' @importFrom glmnet glmnet
#' @importFrom reshape2 melt
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
	fit.lasso <- glmnet(x, y, standardize=T)

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
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=T), stringsAsFactors=FALSE))
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
      df <- as.character(as.data.frame(matrix(unlist(MCR[[i]]), nrow=1, byrow=T), stringsAsFactors=FALSE))
       
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
	fit.lasso <- glmnet(x, y, standardize=T)
	  

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
  CellMorphologyEnrichment <<- as.data.frame(agregatation)
  length(CellMorphologyEnrichment$feature)
  
  #-- Save data --
  return(CMP_subset)
  return(TP_subset)
  return(CellMorphologyEnrichment)

  print("You can find the results in the 'CellMorphologyEnrichment' R objects.")

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