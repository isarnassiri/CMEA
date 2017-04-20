#' @export
#' @import rgl
#' @importFrom clusterSim data.Normalization
#' @import netbenchmark
#' @import gridExtra
#' @importFrom data.table setDT
#' @import glmnet
#' @import vegan
#' @importFrom qgraph centralityTable
#' @importFrom igraph graph.data.frame
#' @importFrom PANR assoScore
#'
#'@title{Mapping query transcriptomic profile against the reference repository}
#'@description{We map the query profile (Q) of fold change of gene expression versus reference repository of transcriptomic data (T) to detect similarities among these profiles (s).}
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'input="BRD-K37798499"
#'Mapping() 
#'}
#'@export

Mapping <- NULL
Mapping <- function() 
{

  #Load data
  data(Transcriptomic_Profile)
  data(Cell_Morphology_Profile)
  
  L1000_TP_profiles <- Transcriptomic_Profile
  L1000_MP_profiles <- Cell_Morphology_Profile
  
  i0 <- which(rownames(L1000_MP_profiles) %in% input)
   
  x_new <- as.data.frame(L1000_TP_profiles[i0,])
  colnames(x_new) <- rownames(L1000_TP_profiles)[i0]
  dim(x_new)
  
  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)
  
  performance <- data.frame()
  TP = TN = FP = FN = 0
  
  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a = query_binary[i]  #my standard
      b = repositoyr_binary[j,i] 
      
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
    TP = TN = FP = FN = 0
  }
  
  colnames(performance) <- c("TP", "TN", "FP", "FN","MCC")
  rownames(performance) <- rownames(repositoyr_binary)
    
  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)
  
  CMP_subset <- L1000_MP_profiles[selected_drugs,]
  TP_subset <- L1000_TP_profiles[selected_drugs,]

  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/CMP_subset.txt", sep = "")
  
  write.table(
    CMP_subset, Destiny_Folder, sep = "\t", row.names = TRUE, quote = FALSE
  )
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/TP_subset.txt", sep = "")
  
  write.table(
    TP_subset, Destiny_Folder, sep = "\t", row.names = TRUE, quote = FALSE
  )
  
  print("You can find the results at: ")
  system.file(package="CMEA")
   
}

#'@export
#'@title{Cell morphology enrichment analysis}
#'@description{We use a stepwise variable selection approach, including combination of least absolute shrinkage and selection operator (LASSO) with cross-validation to tune parameter for cell morphology enrichment analysis. We consider all transcriptomic profiles in the reference repository as inputs of LASSO to select a subset of landmark genes (v) that best describe an indicated profile of cell morphology feature. }
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'input="BRD-K37798499"
#'number_of_features = 20
#'Mapping()
#'Cell_Morphology_Enrichment_Analysis() 
#'}
#'@export

Cell_Morphology_Enrichment_Analysis <- NULL
Cell_Morphology_Enrichment_Analysis <- function()
{
    
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/TP_subset.txt", sep = "")
  TP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/CMP_subset.txt", sep = "")
  CMP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  MCR <- list();
 
  for(i in 1:number_of_features)
  {
    x_new2 <- as.data.frame(CMP_subset[,i])
    colnames(x_new2) <- colnames(CMP_subset)[i]
    
    x <- data.matrix(TP_subset)                #predictors
    y <- as.numeric(unlist(x_new2))            #response
    
    set.seed(1)
    fit.lasso = glmnet(x,y, standardize=TRUE)  
    cv.lasso=cv.glmnet(x,y)
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
  df4 <- (df4[-dim(df4)[1],])  
 
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Temp_file.txt", sep = "")  
  write.table(
    df4, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )

  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Temp_file.txt", sep = "")
  datCM2 <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))   
  Names <- as.character(agregatation$feature)
  
  setDT(agregatation)[, id := .GRP, by = name]    
  agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]  
  agregatation$id <- sprintf("Cluster_%d", agregatation$id)  
  agregatation <- as.data.frame(agregatation)
  length(agregatation$feature)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Cell_Morphology_Enrichment_Analysis_Results.txt", sep = "")
  
  write.table(
    agregatation, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  print("You can find the results at: ")
  system.file(package="CMEA")
     
 }
  
#'@export
#'@title{Ranking of cell morphological phenotypes based on the Strength Centrality Score (SCS).}
#'@description{We consider the connectivity and centrality of cell morphological features in Cosine similarity network of image-based cell morphological profile to rank them based on the Strength Centrality Score (SCS).}
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'TOP=20
#'Ranking_Cell_Morphological_features() 
#'}
#'@export

Ranking_Cell_Morphological_features <- function()
{
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/TP_subset.txt", sep = "")
  TP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/CMP_subset.txt", sep = "")
  CMP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Temp_file.txt", sep = "")
  datCM2 <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))   
  Names <- as.character(agregatation$feature)
  
  setDT(agregatation)[, id := .GRP, by = name]    
  agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]   #sort based on the indices
  agregatation$id <- sprintf("Cluster_%d", agregatation$id) 			#add term of cluster at the beginning of each index
  agregatation <- as.data.frame(agregatation)
  length(agregatation$feature)
  
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
  Long <- Long[!is.na(Long$feature),]
  #g <- ggplot(Long, aes(x = Strength_centrality, y = feature, group = type))
  #g + xlab("") + ylab("") + geom_point() +  geom_path() + labs(list(title = "Top enriched cell morphological phenotypes", x = "Single-cell morphological phenotype score (Strength)", y = "Single-cell morphological phenotype")) + theme_grey(base_size = 10)
 
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Ranking_Cell_Morphological_features.txt", sep = "")
  
  write.table(
    Long, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )
  
  print("You can find the results at: ")
  system.file(package="CMEA")
  
 }
 
#'@export
#'@title{Cross-tabulation of landmark genes, and single-cell morphological features}
#'@description{We present the results of cell morphology enrichment analysis as cross-tabulation of landmark genes, and single-cell morphological features including the direction of effects (up or down regulation).}
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'TOP=10
#'input="BRD-K37798499"
#'crosstabulation() 
#'}
#'@export 
 
crosstabulation <- function()
{
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/TP_subset.txt", sep = "")
  TP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/CMP_subset.txt", sep = "")
  CMP_subset <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Temp_file.txt", sep = "")
  datCM2 <-read.table(Destiny_Folder, sep="\t", header = TRUE)
  
  agregatation <- as.data.frame(aggregate(name ~ feature, data = datCM2, toString))  #aggregation of genes to each CM term
  Names <- as.character(agregatation$feature)
  
  setDT(agregatation)[, id := .GRP, by = name]   						#Package 'data.table' - CRAN, it indexes the gene sets (add new column of indices)
  agregatation <- agregatation[order(agregatation$id, decreasing=FALSE),]   #Sort based on the indices
  agregatation$id <- sprintf("Cluster_%d", agregatation$id)             #Add term of cluster at the beginning of each index
  agregatation <- as.data.frame(agregatation)
  length(agregatation$feature)
  
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
  
  df5 <- datCM2[which(as.character(datCM2$feature) %in% as.character(Long$feature)),]
  df5 <- as.data.frame(df5)
  length(unique(as.character(df5[,2])))
  
  query <- TP_subset[which(rownames(TP_subset) %in% input),]
  
  for(i in 1:dim(df5)[1])
  {
    df5[i,3] <- query[, which(colnames(query) %in% df5[i,1])]
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
  
  df6 <-data.frame()
  
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
  
  MA_f <- MA[1:length(unique(df6[,1])),-(1:length(unique(df6[,1])))]  

  Destiny_Folder <- system.file(package = "CMEA")
  Destiny_Folder = paste(Destiny_Folder, "/Crosstab_table.txt", sep = "")
  
  write.table(
    MA_f, Destiny_Folder, sep = "\t", row.names = TRUE, quote = TRUE
  )
  
  print("You can find the results at: ")
  system.file(package="CMEA")

 }
 
#'@export
#'@title{Modeling of cell morphological features based on the transcriptomic profile.}
#'@description{We leverage the significant cross correlation between the cell morphological feature and related transcriptomic profiles to predict previously unrecognized cell morphological states for transcriptomic of experimental perturbation of interest.
#'You should use the whole of transcriptomic and cell morphological profiles as input of this function to get real results.}
#'@author{Isar Nassiri, Matthew McCall}
#'@examples{
#'data(Transcriptomic_Profile)
#'data(Cell_Morphology_Profile)
#'Number_features=5  #Number of cell morphological phenotype
#'Number_profiles=5  #An integer, number of profiles (test sets)
#'Modeling_morphological_features() 
#'}
#'@export 

Modeling_morphological_features <- NULL
Modeling_morphological_features <- function()
{
#Load data
data(Transcriptomic_Profile)
data(Cell_Morphology_Profile)

L1000_TP_profiles <- Transcriptomic_Profile
L1000_MP_profiles <- Cell_Morphology_Profile

Number_profiles = 5
Number_features = 5
rand <- sample(1:dim(L1000_MP_profiles)[2], Number_features, replace = FALSE)

a1 <- data.frame()

for(i0 in 1:Number_profiles)
{
  x_new <- as.data.frame(L1000_TP_profiles[i0,])
  colnames(x_new) <- rownames(L1000_TP_profiles)[i0]
  
  query_binary <- ifelse(x_new > 0, 1, 0)
  repositoyr_binary <- ifelse(L1000_TP_profiles > 0, 1, 0)
  
  performance <- data.frame()
  TP = TN = FP = FN = 0
  
  for(j in 1:dim(repositoyr_binary)[1])
  {
    for(i in 1:dim(repositoyr_binary)[2])
    {
      a = query_binary[i]   
      b = repositoyr_binary[j,i] 
      
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
    TP = TN = FP = FN = 0
  }
  
  colnames(performance) <- c("TP", "TN", "FP", "FN","MCC")
  rownames(performance) <- rownames(repositoyr_binary)
  
  selected_drugs <- which(performance[,5] > 0.1)
  length(selected_drugs)
  
  CMP_subset <- L1000_MP_profiles[selected_drugs,] 
  TP_subset <- L1000_TP_profiles[selected_drugs,]
  
  CMP_subsetn <- data.Normalization(CMP_subset,type="n4");
  TP_subsetn <- data.Normalization(TP_subset,type="n4");
  CMP_subset0 <- data.matrix(CMP_subsetn)
  TP_subset0 <- data.matrix(TP_subsetn)

  pred1 <- data.frame()
  
  y0 <- 1
  
  for(j in 1:Number_features)
  {
    
    i <- which(rownames(TP_subset0) %in% rownames(L1000_TP_profiles)[i0])
    
    x = TP_subset0[-i,]
    y = CMP_subset0[-i,rand[j]]
    
    lasso.2 <- glmnet(x, y, standardize=TRUE)
    
    #-- extract significant coefficients -- 
    train=sample(seq(dim(x)[1]),(dim(x)[1]/2),replace=FALSE)
    lasso.tr=glmnet(x[train,],y[train])
    pred=predict(lasso.tr,x[-train,])
    rmse= sqrt(apply((y[-train]-pred)^2,2,mean))
    lam.best=lasso.tr$lambda[order(rmse)[1]]
    lam.best
    
    Lasso_coefficient <- coef(lasso.2, s=lam.best)
    inds<-which(Lasso_coefficient[,1]!=0)
    variables<-row.names(Lasso_coefficient)[inds]
    variables<-variables[!(variables %in% '(Intercept)')];
    
    length(which(0 != Lasso_coefficient[,1]))
    results <- as.data.frame(Lasso_coefficient[which(0 != Lasso_coefficient[,1]),1])
    colnames(results) <- "coef"
    results <- as.data.frame(results[-1,,FALSE])  # FALSE is about inactivation of drop paremters
    dim(results)[1]
    
    if(1<dim(results)[1])
    {			
      pred <- predict(lasso.2, newx=TP_subset0[i,,drop = FALSE], s=lam.best)
      
      pred1[y0,1] <- as.numeric(pred)
      pred1[y0,2] <- CMP_subset0[i, rand[j]]
      pred1[y0,3] <- colnames(CMP_subset0)[rand[j]]
      y0 = y0 + 1
    }
  }
  
  if(0<dim(pred1)[1])
  {
    colnames(pred1) <- c("Prediction", "Experiment", "CM")
    pred1 <- pred1[!is.na(pred1$Prediction),]
    a1 <- rbind(a1, pred1)
  }
  
}

name <- unique(a1$CM)

a2 <- data.frame()

par(mfrow=c(2,5))

for(i in 1:length(name))
{
  Experiment = a1$Experiment[which(a1$CM %in% name[i])]
  Prediction = a1$Prediction[which(a1$CM %in% name[i])]
  
  a2[i,1] <- name[i]
  
  a2[i,2] <- (cor(a1$Experiment[which(a1$CM %in% name[i])], 
                  a1$Prediction[which(a1$CM %in% name[i])]))^2
  
}

	colnames(a2) <- c("Cell morphology", "Exp. and Pred. Correlations (R^2)")

	# -- Save data --
	Destiny_Folder <- system.file(package = "CMEA")
	Destiny_Folder = paste(Destiny_Folder, "/Modeling_correlation.txt", sep = "")

	write.table(
	  a2, Destiny_Folder, sep = "\t", row.names = TRUE, quote = TRUE
	)

	Destiny_Folder <- system.file(package = "CMEA")
	Destiny_Folder = paste(Destiny_Folder, "/Modeling_details.txt", sep = "")

	write.table(
	  a1, Destiny_Folder, sep = "\t", row.names = TRUE, quote = TRUE
	)

	print("You can find the results at: ")
	system.file(package="CMEA") 

 }