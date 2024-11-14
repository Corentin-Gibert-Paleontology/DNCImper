#' DNCI Function : Dispersal-Niche Continuum Index computation for overall dissimilarity analysis of 3 groups or more
#'
#' Quantitative identification of the main assembly process driving the overall dissimilarity contribution distribution (i.e. the empirical SIMPER profile structure).
#' This function is based on PerSIMPER_overall function and its E index return(). -Under development-
#' The three distributions of E index (corresponding to the three hypothesis: niche, dispersal, niche+dispersal) are used to compute the DNCI index.
#' If DNCI is significantly < 0 : dispersal || DNCI significantly > 0 : niche || DNCI +- CI ~~ 0 : dispersal+niche
#' -under development-
#' @param x Sample/Taxa matrix with sample in row and taxa in column
#' @param grouping Grouping vector, ex : c(1,1,1,1,2,2,2,2,2,3,3,3) : 3 groups or more
#' @param id Name of the dataset, default = "no_name"
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataTYPE Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @param parallelComputing Run PerSIMPER on half of the available cores/nodes
#' @examples A <- DNCImper:::DNCI.ses_overall(DNCImper::Matrix_4groups, DNCImper::Group4)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 3 groups or more.
#' @examples #
#' @examples B <- DNCImper:::DNCI.ses_overall(DNCImper::Matrix_4groups, DNCImper::Group4, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
#' @importFrom graphics legend lines title
#' @importFrom stats median quantile sd
#' @importFrom utils combn
#'
#'





################ UNDER DEVELOPMENT ######################
## To contact me : corentingibert@gmail.com (feel free)
################ EN DEVELOPPEMENT #######################
# This function calls/needs PerSIMPER_Overall function #


################ UNDER DEVELOPMENT ##################
####### DNCI.ses function modified to be used #######
####### on multiples groups matrix (more than 2) ####
###### to compute OVERALL PerSIMPER/DNCI analysis ###
###### OVERALL = contribution of taxa to the ########
###### overall between-group dissimilarity ##########
###### and not to individual pair of groups #########
#####################################################

##### Same arguments and results as in DNCI_ses.R

DNCI.ses_overall <- function(x, grouping, id = "no_name", Nperm = 1000, count = TRUE, dataTYPE = "prab",
                                     plotSIMPER = TRUE, parallelComputing = FALSE) { #this calculates the metric using PERSIMPER - now the output included DELTAd-n, sd of DELTA.d-n and confidence interval
  groups <- sort(unique(grouping))
  #results = PerSIMPER(x, grouping, log = TRUE, count = count)
  results = DNCImper:::PerSIMPER_overall(x, grouping, count = count, Nperm = Nperm, dataTYPE = dataTYPE, plotSIMPER = plotSIMPER, parallelComputing = parallelComputing)
  E = results[["EcartCarreLog"]]

  #first calculate SES.d and SES.n based on E values from PERSIMPER
  SES.d = (E$Orange-mean(E$Blue))/sd(E$Blue)
  SES.n = (E$Green-mean(E$Blue))/sd(E$Blue)

  #then calculate DELTA.d-n
  DELTA.dn = mean(SES.d)-mean(SES.n)
  #then calculate sd related to DELTA.d-n
  S.DELTA.dn = sqrt(sd(SES.d)^2+sd(SES.n)^2)
  #get the confidence interval based on S.DELTA.dn
  CI.DELTA.dn = 2*S.DELTA.dn
  # ---> then you have both the DELTA.d-n and the 95% confidence interval wich is 2 * S.DELTA.d-n
  metric = data.frame(id = id, group1= groups[1], group2 = groups[2],DELTA.dn, CI.DELTA.dn, S.DELTA.dn)
  #,id dans le function() a ?t? supprim? aussi

  return(metric)
}
