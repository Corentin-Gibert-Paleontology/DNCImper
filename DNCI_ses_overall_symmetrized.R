
#' DNCI Function : Dispersal-Niche Continuum Index computation for overall dissimilarity analysis of 3 groups or more with randomly even groups
#'
#' Quantitative identification of the main assembly process driving the overall dissimilarity contribution distribution (i.e. the empirical SIMPER profile structure).
#' This function is based on DNCI.ses_overall() PerSIMPER_overall function and its E index return(). -Under development-
#' -under development-
#' Group are made even by subsampling largest group to the size of the smallest ! CAUTION ! Do rerun for robust results
#' -under development-
#' @param Mat Sample/Taxa matrix with sample in row and taxa in column
#' @param Group Grouping vector, ex : c(1,1,1,1,2,2,2,2,2,3,3,3) : 3 groups or more
#' @param id Name of the dataset, default = "no_name"
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataTYPE Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @param NbrReRun Number of iteration to obtain mean DNCI_overall values with even groups
#' @param parallelComputing Run PerSIMPER_overall on half of the available cores/nodes
#' @examples A <- DNCImper:::DNCI.ses_overall_symmetrized(DNCImper::Matrix_4groups, DNCImper::Group4, NbrReRun = 10)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 3 groups or more; 10 reruns (default NbrReRun is 100).
#' @examples #
#' @examples B <- DNCImper:::DNCI.ses_overall_symmetrized(DNCImper::Matrix_4groups, DNCImper::Group4, NbrReRun = 10, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations and 10 reruns (default NbrReRun is 100), with no countdown and no plots
#' @importFrom graphics legend lines title
#' @importFrom stats median quantile sd
#' @importFrom utils combn
#'
#'

################ UNDER DEVELOPMENT ######################
## To contact me : corentingibert@gmail.com (feel free)
################ EN DEVELOPPEMENT #######################


# This function calls/needs DNCI_ses_overall function ###
# This function calls/needs PerSIMPER_Overall function ##

######################### CAUTION #######################
#### This function is a wraper for DNCI_ses_overall to ##
### compute robust overall DNCI values by making all ####
### groups even. ########################################
##### By making all groups even, this function can ######
##### change massively the original DNCI_overall ########
###### IF (!) one group is way smaller than the others ##
############## Consequently you should (A) run a lot ####
########### of iterations or (B) remove the smallest ####
######## group ##########################################
#########################################################

### Arguments are similar to DNCI_ses_overall.R
##### BUT (!) you can choose the number of iteration to
######## compute mean robust overall DNCI value : NbrRerun

DNCI.ses_overall_symmetrized <- function(Mat, Group, id = "no_name", NbrReRun = 100, dataTYPE = "prab",
                                                 Nperm = 100, plotSIMPER = TRUE, count = FALSE, parallelComputing = FALSE)
{
  warning("This function is a wrapper of PerSIMPER_overall and DNCI_overall functions")
  warning("It is still under developement")
  warning("Please use NbrReRun argument to obtain mean/median results based on multiple resampling")
  warning("Keep in mind that large group will be resampled to the size of the smallest one")
  warning("If a large gap exist between groups size, please consider large number of replicate (NbrReRun) and/or exclude the smallest group")
  Number_ofpairs <- combn(1:length(unique(Group)), 2)
  Unik_group <- unique(Group)

  List_result <- c()

  for(x in 1:NbrReRun)
  {
    SampleGroup <- c()
    for(y in 1:length(Unik_group))
    {
      SampleGroup <- c(SampleGroup, sample(which(Group == y), min(table(Group))))
    }
    Mat_Sampled <- Mat[SampleGroup,]
    Group_Sampled <- Group[SampleGroup]
    if(length(which(apply(Mat_Sampled, 2, sum) == 0)) != 0)
    {
      Mat_Sampled <- Mat_Sampled[,-which(apply(Mat_Sampled, 2, sum) == 0)]
    }

    if(length(Number_ofpairs[1,]) > 2){

      Analyse_Overall <- DNCImper:::DNCI.ses_overall(Mat_Sampled, Group_Sampled, id = id, Nperm = Nperm, plotSIMPER = plotSIMPER,
                                                 count = count, dataTYPE = dataTYPE, parallelComputing = parallelComputing)

      if(x == 1)
      {
        List_result <- Analyse_Overall
      }
      if(x != 1)
      {
        List_result <- rbind(List_result, Analyse_Overall)
      }
    }
  }
  return(List_result)
}
