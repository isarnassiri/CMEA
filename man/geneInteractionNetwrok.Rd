% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CMEA.R
\name{geneInteractionNetwrok}
\alias{geneInteractionNetwrok}
\title{Gene regulatory network of cell morphological phenotypes}
\usage{
geneInteractionNetwrok(number_of_features, lift, confidence,
  Query_Transcriptomic_Profile, Transcriptomic_Profile, Cell_Morphology_Profile)
}
\arguments{
\item{number_of_features}{The number of top cell morphological features, ranked based on the Strength Centrality Score (SCS) for enrichment analysis (e.g. 20).
This parameter specifies the number of first top cell morphological features which are used for enrichment sets of landmark genes.}

\item{lift}{We use an association mining model to detect the associations between the genes and to infer a gene regulatory network from phenotypic experimental data. In order to select significant interactions between any two genes from
the complete digraph of interactions between genes, we use minimum thresholds on lift and confidence measures. Lift measures the frequency of an indicated gene in the cell morphological gene sets.}

\item{confidence}{We use an association mining model to detect the associations between the genes and to infer a gene regulatory network from phenotypic experimental data. In order to select significant interactions between any two genes from
the complete digraph of interactions between genes, we use minimum thresholds on lift and confidence measures. Confidence shows how often a given association
rule between two genes has been found in the dataset.}

\item{Query_Transcriptomic_Profile}{A data frame including expression level of 978 land mark genes in response to treatment with an indicated drugs/small compound molecule (row drug/small molecule compound ID, and column land mark gene symbols)}

\item{Transcriptomic_Profile}{A data frame including expression level of 978 land mark genes in response to treatment with 162 drugs/small compound molecules (row drug/small molecule compound ID, and column land mark gene symbols)}

\item{Cell_Morphology_Profile}{A data frame including profiles of 812 cell morphological features of 162 drugs/small compound molecules (row drug/small molecule compound ID, and column cell morphological features)}
}
\value{
1. A matrix of gene-gene interaction network, including two columns under title of Source and Target. Source molecule in network interaction with Target molecule based on biological rationale.
}
\description{
We use the mining association model to detect the associations between the landmark genes, and inference of gene regulatory network of cell morphological phenotypes.
}
\examples{
data(Transcriptomic_Profile)
data(Cell_Morphology_Profile)
data(Query_Transcriptomic_Profile)
geneInteractionNetwrok(10, 0.1, 0.6, Query_Transcriptomic_Profile, Transcriptomic_Profile, Cell_Morphology_Profile)
}
\author{
{Isar Nassiri, Matthew McCall}
}
