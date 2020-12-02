#' DNCI multigroup Function : Pairwise Dispersal-Niche Continuum Index on 2 groups or more
#'
#' This function is a wraper for DNCI.ses() function. Can be used with 2 group or more.
#'
#' Quantitative identification of the main assembly process. Pairwise analysis.
#' This function is based on DNCI_ses (for 2 group analysis) and PerSIMPER function E index return().
#' The three distributions of E index (corresponding to the three hypothesis: niche, dispersal, niche+dispersal) are used to compute the DNCI index.
#' If DNCI is significantly < 0 : dispersal || DNCI significantly > 0 : niche || DNCI +- CI ~~ 0 : dispersal+niche
#' See Vilmi et al. 2021 Ecography for DNCI computation and more information on process identification as well as example with Chinese and US river communities.
#' More information in code and comments inside function file.
#' @param x Sample/Taxa matrix with sample in row and taxa in column
#' @param grouping Grouping vector, ex : c(1,1,1,1,2,2,2,2,2) : 2 groups or more ex2 : c(1,1,1,1,1,2,2,2,2,3,3,3,3)
#' @param id Name of the dataset, default = "no_name"
#' @param symmetrize By random sampling of the largest groups, analyzed pairs of group are made even. Strong heterogeneity in group length can impact DNCI values. More information in function file and Vilmi et al. 2021.
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataType Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @examples A <- DNCI_multigroup(Matrix, Group)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 2 groups ONLY
#' @examples #
#' @examples B <- DNCI_multigroup(Matrix, Group, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
#'
#'



#############################################################
######### DNCI Function for 2 groups AND more  ##############
################# ANALYSIS ON PAIRS #########################
###### This function calls PerSIMPER and DNCI functions #####

## x = matrix with taxa in columns and localities in rows
## grouping = grouping vector for localities e.g. group <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2)
## grouping need to have the same length as X number of rows
## id = Name of your dataset
## Nperm and count arguments are for PerSIMPER calling, same argument as PerSIMPER fun()

## SYMMETRIZE : [IMPORTANT] : this argument make group even by subsampling the largest
##                            group to reduce it (or them) to the smallest group size
##              [IMPORTANT] : Repeat computation X (e.g. 1000) times to obtain mean values
##                            Effect can be strong if groups are strongly uneven

DNCI_multigroup <- function(x, grouping,id = "no_name", Nperm = 1000, count = TRUE, symmetrize = FALSE, plotSIMPER = TRUE) {
  group.combinations <- combn(unique(sort(grouping)),2)

  ddelta <- NULL

  if(NCOL(group.combinations) == 1) #If only 2 groups are compared
  {
    ddelta <- DNCI.ses(x , grouping,id=id, Nperm = Nperm, count = count, plotSIMPER = plotSIMPER)
  }
  if(NCOL(group.combinations) > 1)
  {

    for(i in 1:NCOL(group.combinations)) {
      splitx <- split(x,grouping)

      #Ici symmetrize:
      if(symmetrize == TRUE)
      {
        Add <- which(c(NROW(splitx[[group.combinations[1,i]]]),
                       NROW(splitx[[group.combinations[2,i]]])) == max(c(NROW(splitx[[group.combinations[1,i]]]),
                                                                         NROW(splitx[[group.combinations[2,i]]]))))
        if(Add == 1)
        {

          sampled_lines <- sample(1:length(splitx[[group.combinations[1,i]]][,1]),
                                  length(splitx[[group.combinations[2,i]]][,1]))
          splitx[[group.combinations[1,i]]] <- splitx[[group.combinations[1,i]]][sampled_lines,]
        }

        if(Add == 2)
        {

          sampled_lines <- sample(1:length(splitx[[group.combinations[2,i]]][,1]),
                                  length(splitx[[group.combinations[1,i]]][,1]))
          splitx[[group.combinations[2,i]]] <- splitx[[group.combinations[2,i]]][sampled_lines,]
        }

      }

      paired.x <- rbind(splitx[[group.combinations[1,i]]],
                        splitx[[group.combinations[2,i]]])

      # remove empty species
      ifzero <- which(apply(paired.x, 2, sum) == 0)
      if(length(ifzero > 0)){
        paired.x <- paired.x[,-which(colSums(paired.x)==0)]}
      if(length(which(rowSums(paired.x) == 0)) != 0){stop("ERROR : A row/sample is empty")}
      group.pair <- c(rep(group.combinations[1,i], NROW(splitx[[group.combinations[1,i]]])),
                      rep(group.combinations[2,i], NROW(splitx[[group.combinations[2,i]]])))
      ddelta <- rbind(ddelta, DNCI.ses(x=paired.x,grouping=group.pair,id=id, Nperm = Nperm, count = count, plotSIMPER = plotSIMPER)) #here is the part that calculates the index based on PERSIMPER
    }
  }
  return(ddelta)
}

### return : similar to DNCI results with the exception that returned variable is a multiples rows dataframe().
