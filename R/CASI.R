#' Canonical Analysis of Set Interactions
#'
#' Assuming a case-control type study design, two sets of features are measured in both cases and controls.
#' The objective of CASI hypothesis testing is to evaluate evidence of statistical interactions between these 
#' sets in relation to outcome status. 
#'
#' The two feature sets could be sets of SNPs underlying two genes, a set of environmental exposures and a set of 
#' SNPs, a set of SNPs and a set of DNA methylation probes, etc. The null hypothesis is that there are no two linear 
#' combinations of the two sets whose correlation differs between cases and controls. Such a difference in correlation
#' would imply a statistical interaction. The distribution of the CASI statistic is unknown, so rather than provide a 
#' p-value, the CASI function generates the test statistic under the observed and null hypothesis conditions. The null 
#' conditions are imposed by randomly permuting case-control labels and repeating the analysis n.perm times. It is
#' anticipated by the developers that an FDR approach will be applied to CASI function results generated from 
#' applications to multiple, perhaps thousands of pairs of features. However, a permutation-based p-value could be
#' computed if the number of permutations is sufficiently large. 
#'
#' @param G1_cs variable set 1 measured in cases. Variables are organized in columns and individuals in rows. 
#' Therefore, the number of rows equals the number of cases.
#' @param G2_cs variable set 2 measured in cases. Variables are organized in columns and individuals in rows. 
#' Therefore, the number of rows equals the number of cases.
#' @param G1_cn variable set 1 measured in controls. Variables are organized in columns and individuals in rows. 
#' Therefore, the number of rows equals the number of controls.
#' @param G2_cn variable set 2 measured in controls. Variables are organized in columns and individuals in rows. 
#' Therefore, the number of rows equals the number of controls.
#' @param lst.perm a list with length equal to the number of permutations to be conducted. Elements are vectors of  
#'  intergers, each of length equal to the total sample size (n = n.cases + n.controls). Each vector contains indices
#' specifying a random permutation for case-control labels, where the observed order is represented by vector 1:n. 
#' The same lst.perm should be used for multiple tests to capture possible dependencies among tests for use in fdrci. 
#' For example: lst.perm[[1]] <- sample(1:n), lst.perm[[2]] <- sample(1:n), ...
#' @param fs typically the larger the CASI statistic, the more extreme. However, if fs is set to TRUE the CASI statistic 
#' is flipped so that the smaller it is, the more extreme it is, like a p-value. This facilitates the use of the fdrci R 
#' package for computing FDR estimates and corresponding confidence intervals.
#'
#' @return A vector of length n.perm + 1, where the first element of the vector is the CASI statistic from the observed
#' data and the following n.perm elements are CASI statistics computed using the permuted data.
#' 
#' @author Joshua Millstein, \email{joshua.millstein@usc.edu}, Vladimir Kogan
#' @references Vladimir Kogan and Joshua Millstein. 2018. Genetic-Epigenetic Interactions in Asthma Revealed by a 
#' Genome-Wide Gene-Centric Search. Human Heredity (in review)
#' @importFrom PMA CCA
#' @importFrom stats sd
#' @keywords GxG interactions, GxE interactions, canonical correlation analysis, case-only interaction analysis
#'
#' @examples
#' n.case = 100
#' n.control = 100
#' m.set1 = 10
#' m.set2 = 10
#' n.effects = 5
#' beta.vec = runif(n.effects, 1, 2)
#' set1.nms = paste("V1.", 1:m.set1, sep="")
#' set2.nms = paste("V2.", 1:m.set1, sep="")
#' nms = c("case", set1.nms, set2.nms)
#' mydat = as.data.frame(matrix(NA, nrow=0, ncol=length(nms)))
#' names(mydat) = nms
#' logit = function(p) log(p / (1-p))
#' logistic = function(a) exp(a) / (exp(a)+1)
#' risk.base = .05
#' 
#' rowno = 0
#' ind.case = 0
#' ind.control = 0
#' while(rowno < (n.case+n.control)){
#' 	vec1 = rbinom(m.set1, 2, .2)
#' 	vec2 = rbinom(m.set2, 2, .2)
#' 	lc = logit(risk.base)
#' 	for(i in 1:n.effects) lc = lc + beta.vec[i]*vec1[i]*vec2[i]
#' 	myrisk = logistic(lc)
#' 	if( runif(1) < myrisk ){ 
#' 		if(ind.case < n.case){ 
#' 			ind.case = ind.case + 1
#' 			rowno = ind.case + ind.control
#' 			mydat[ rowno,"case"] = 1
#' 			mydat[ rowno, c(set1.nms, set2.nms) ] = c(vec1, vec2)
#' 		}
#' 	} else { 
#' 		if(ind.control < n.control){ 
#' 			ind.control = ind.control + 1
#' 			rowno = ind.case + ind.control
#' 			mydat[ rowno,"case"] = 0
#' 			mydat[ rowno, c(set1.nms, set2.nms) ] = c(vec1, vec2)
#' 		}
#' 	}
#' } # end while()
#' 
#' case.status = mydat[, "case"]
#' G1_cs = mydat[ is.element(case.status, 1), set1.nms ]
#' G2_cs = mydat[ is.element(case.status, 1), set2.nms ]
#' G1_cn = mydat[ is.element(case.status, 0), set1.nms ]
#' G2_cn = mydat[ is.element(case.status, 0), set2.nms ]
#' n.perm = 10
#' lst.perm = vector('list', n.perm)
#' for(p in 1:n.perm) lst.perm[[p]] = sample(1:(n.case+n.control))
#' CASI(G1_cs, G2_cs, G1_cn, G2_cn, lst.perm)
#'
#' @export
CASI <- function(G1_cs, G2_cs, G1_cn, G2_cn, lst.perm, fs=FALSE){

   diff_fish_obs<-try(diff_fish(G1_cs, G2_cs, G1_cn, G2_cn)[1], silent = T)

   if(class(diff_fish_obs) %in% "try-error"){
      diff_fish_obs <- try(diff_fish_SCCA(G1_cs, G2_cs, G1_cn, G2_cn)[1], silent = T)
   }  
   G1_names <- names(G1_cs)
   G2_names <- names(G2_cs)
   G1_2_cs <- cbind(G1_cs, G2_cs)
   G1_2_cs$CaseStatus <- "case"
   G1_2_cn <- cbind(G1_cn, G2_cn)
   G1_2_cn$CaseStatus <- "control"
   mydat <- rbind(G1_2_cs, G1_2_cn)

   # Permute case-control labels and recalculate CASI length(n.perm) times
   CASI_perm <- vector()
   for(i in 1:length(lst.perm)){
      mydat$CaseStatus <- mydat$CaseStatus[lst.perm[[i]]]
 
      G1_cs <- subset(mydat, CaseStatus=="case")[,G1_names] 
      G2_cs <- subset(mydat, CaseStatus=="case")[,G2_names]  
  
      G1_cn <- subset(mydat, CaseStatus=="control")[,G1_names]  
      G2_cn <- subset(mydat, CaseStatus=="control")[,G2_names] 
  
      # Calculate canonical correlation and fisher transformed difference
      diff_fish_perm <- try(diff_fish(G1_cs, G2_cs, G1_cn, G2_cn), silent = T)
  
      if(class(diff_fish_perm) %in% "try-error"){
          diff_fish_perm <- try(diff_fish_SCCA(G1_cs, G2_cs, G1_cn, G2_cn), silent = T)
      }  
      CASI_perm[i] <- diff_fish_perm[1]
   }

   #  standardize observed
   diff_obs_std <- (diff_fish_obs-mean(CASI_perm))/sd(CASI_perm)

   # standardize permuted
   diff_perm_std <- (CASI_perm-mean(CASI_perm))/sd(CASI_perm)

   CASI_stats <- c(diff_obs_std, diff_perm_std)
   if(fs) CASI_stats <- 10^-abs(CASI_stats) # flip to make smaller-is-better for fdrci functions
   names(CASI_stats) <- c("Observed", paste("Perm", 1:length(lst.perm)))

   return(CASI_stats)
 
} # end CASI()

# fisher's z transformation
fisherz <- function (rho) {
    0.5 * log((1 + rho)/(1 - rho))
}

# function for calculating fisher transformed difference using the default canonical correlation function    
diff_fish <- function(G1_cs, G2_cs, G1_cn, G2_cn){  
  
   # canonical correlation for cases 
   cc_cs <- cancor(G1_cs, G2_cs)
  
   # largest canonical correlation
   cs_cc <- cc_cs$cor[1]

   # coefficients for Variable Set 1 
   G1_cs_coefs <- cc_cs$xcoef[,1] 
   # coefficients for Variable Set 2 
   G2_cs_coefs <- cc_cs$ycoef[,1]
  
   # canonical variate for controls for variable set 1 
   G1_cn_v <- as.matrix(G1_cn[,attributes(G1_cs_coefs)$names]) %*% G1_cs_coefs

   # canonical variate for controls for variable set 2
   G2_cn_v <- as.matrix(G2_cn[,attributes(G2_cs_coefs)$names]) %*% G2_cs_coefs
  
   # correlation between control canonical variates
   cn_cc<-cor(G1_cn_v, G2_cn_v)
  
   # difference between fisher transformed case and control canonical variates
   fz_diff = as.numeric(fisherz(cs_cc)-fisherz(cn_cc))
   casi_rslt_vec <- c(fz_diff, cs_cc, cn_cc)
 
   return(casi_rslt_vec)
   
} # end diff_fish()

# function for calculating fisher transformed difference using the Sparse canonical correlation function    
diff_fish_SCCA <- function(G1_cs, G2_cs, G1_cn, G2_cn){
  
  # Canonical Correlation Observed
  cc_cs <- PMA::CCA(G1_cs, G2_cs, standardize=T, K=4)
  cs_cc <- max(cc_cs$cors)
  cc_cs_v1_coefs <- cc_cs$u[,match(cs_cc, cc_cs$cors)] 
  cc_cs_v2_coefs <- cc_cs$v[,match(cs_cc, cc_cs$cors)]
  cn_v1_v <- as.matrix(G1_cn) %*% cc_cs_v1_coefs
  cn_v2_v <- as.matrix(G2_cn) %*% cc_cs_v2_coefs
  cn_cc <- cor(cn_v1_v, cn_v2_v)
  
  fz_diff <- as.numeric(fisherz(cs_cc) - fisherz(cn_cc))
  
  casi_rslt_vec <- c(fz_diff, cs_cc, cn_cc)
  
  return(casi_rslt_vec)
  
} # end diff_fish_SCCA()


