#' @export
#' @importFrom clusterSim data.Normalization
#' @importFrom netbenchmark clr.wrap
#' @importFrom data.table setDT
#' @importFrom glmnet glmnet
#'

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