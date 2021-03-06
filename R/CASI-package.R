#' Canonical Analysis of Set Interactions
#' 
#' Assuming a case-control type study design, two sets of features are measured in both cases and controls.
#' The objective of CASI hypothesis testing is to evaluate evidence of statistical interactions between these 
#' sets in relation to outcome status. 
#' 
#' \tabular{ll}{ Package: \tab CASI\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2018-04-28\cr License: \tab Artistic-2.0\cr LazyLoad: \tab
#' yes\cr } Assuming a case-control type study design, two sets of features are measured in both cases and controls.
#' The objective of CASI hypothesis testing is to evaluate evidence of statistical interactions between these 
#' sets in relation to outcome status. The two feature sets could be sets of SNPs underlying two genes, a set of 
#' environmental exposures and a set of SNPs, a set of SNPs and a set of DNA methylation probes, etc. The null 
#' hypothesis is that there are no two linear combinations of the two sets whose correlation differs between cases 
#' and controls. Such a difference in correlation would imply a statistical interaction. The distribution of the CASI 
#' statistic is unknown, so rather than provide a p-value, the CASI function generates the test statistic under the 
#' observed and null hypothesis conditions. The null conditions are imposed by randomly permuting case-control 
#' labels and repeating the analysis n.perm times. It is anticipated by the developers that an FDR approach will be 
#' applied to CASI function results generated from applications to multiple, perhaps thousands of pairs of features. 
#' However, a permutation-based p-value could be computed if the number of permutations is sufficiently large. 
#' 
#' @name CASI-package
#' @aliases CASI-package CASI
#' @docType package
#' @author Joshua Millstein, Vladimir Kogan
#' 
#' Maintainer: Joshua Millstein <joshua.millstein@usc.edu> Joshua Millstein
#' @references Vladimir Kogan and Joshua Millstein. 2018. Genetic-Epigenetic Interactions in Asthma Revealed by a 
#' Genome-Wide Gene-Centric Search. Human Heredity (in review)
#' @keywords GxG interactions, GxE interactions, canonical correlation analysis, case-only interaction analysis
NULL



