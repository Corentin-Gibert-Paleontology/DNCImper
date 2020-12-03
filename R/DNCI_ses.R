#' DNCI Function : Dispersal-Niche Continuum Index
#'
#' Quantitative identification of the main assembly process.
#' This function is based on PerSIMPER function and its E index return().
#' The three distributions of E index (corresponding to the three hypothesis: niche, dispersal, niche+dispersal) are used to compute the DNCI index.
#' If DNCI is significantly < 0 : dispersal || DNCI significantly > 0 : niche || DNCI +- CI ~~ 0 : dispersal+niche
#' See Vilmi, Gibert et al. 2021 Ecography for DNCI computation and more information on process identification.
#' More information in code and comments inside function file.
#' @param x Sample/Taxa matrix with sample in row and taxa in column
#' @param grouping Grouping vector, ex : c(1,1,1,1,2,2,2,2,2) : 2 groups only !!
#' @param id Name of the dataset, default = "no_name"
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataType Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @examples A <- DNCI.ses(Matrix, Group)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 2 groups ONLY
#' @examples #
#' @examples B <- DNCI.ses(Matrix, Group, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
#'
#'

#############################################################
###### DNCI Function : Dispersal-Niche Continuum Index ######
######### idea by Gibert, Escarguel, Vilmi and Wang #########
#############################################################
###### This function calls/needs PerSIMPER function #########

##### DIFFERENCES WITH PerSIMPER : QUANTITATIVE Results #####

## x = matrix with taxa in columns and localities in rows
## grouping = grouping vector for localities e.g. group <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2)
## grouping need to have the same length as X number of rows
## id = Name of your dataset
## Nperm and count arguments are for PerSIMPER calling, same argument as PerSIMPER fun()

## To contact me : corentingibert@gmail.com or annika.vilmi@gmail.com (feel free)

DNCI.ses <- function(x, grouping,id = "no_name", Nperm = 1000, count = TRUE, plotSIMPER = TRUE) { #this calculates the metric using PERSIMPER - now the output included DELTAd-n, sd of DELTA.d-n and confidence interval
  groups <- sort(unique(grouping))
  stopifnot(length(groups) == 2)
  results = PerSIMPER(x, grouping,  count = count, Nperm = Nperm, plotSIMPER = plotSIMPER)
  E = results[["EcartCarreLog"]]

  #first calculate SES.d and SES.n based on E values from PERSIMPER
  SES.d = (E$Orange-mean(E$Blue))/sd(E$Blue)
  SES.n = (E$Green-mean(E$Blue))/sd(E$Blue)

  #then calculate DNCI
  DNCI = mean(SES.d)-mean(SES.n)
  #then calculate sd related to DNCI
  S.DNCI = sqrt(sd(SES.d)^2+sd(SES.n)^2)
  #get the confidence interval based on S.DNCI
  CI.DNCI = 2*S.DNCI
  # ---> then you have both the DNCI and the 95% confidence interval wich is 2 * S.DNCI
  metric = data.frame(id=id, group1= groups[1], group2 = groups[2], DNCI, CI.DNCI, S.DNCI)


  return(metric)
}

### return : DNCI = Dispersal-Niche continuum index value. Negative values = Dispersal ; Positive values = Niche
###          CI.DNCI = Confidence interval associated with DNCI; mandatory to read DNCI analysis
###          S.DNCI = Variance associated with DNCI
