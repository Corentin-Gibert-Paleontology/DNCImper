
#' DNCI Function : Dispersal-Niche Continuum Index computation for overall dissimilarity analysis of 3 groups or more with randomly even groups
#'
#' Quantitative identification of the main assembly process driving the overall dissimilarity contribution distribution (i.e. the empirical SIMPER profile structure).
#' This function is based on DNCI.ses_overall() PerSIMPER_overall function and its E index return(). -Under development-
#' -under development-
#' Group are made even by subsampling largest group to the size of the smallest ! CAUTION ! Do rerun for robust results
#' -under development-
#' @param x Sample/Taxa matrix with sample in row and taxa in column
#' @param grouping Grouping vector, ex : c(1,1,1,1,2,2,2,2,2,3,3,3) : 3 groups or more
#' @param id Name of the dataset, default = "no_name"
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataType Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @param NbrReRun Number of iteration to obtain mean DNCI_overall values with even groups
#' @examples A <- DNCI.ses_overall(Matrix, Group)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 3 groups or more.
#' @examples #
#' @examples B <- DNCI.ses_overall(Matrix, Group, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
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

DNCI.ses_overall_symmetrized <- function(Mat, Group, id = "no_name", NbrReRun = 100, Nperm = 100, count = FALSE)
{
  Number_ofpairs <- combn(1:length(unique(Group)), 2)
  Unik_group <- unique(Group)

  List_result <- c()

  for(x in 1:NbrReRun)
  {
    SampleGroup <- c()
    for(y in 1:length(Unik_group)) #CHANGER PAR LA
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
      #Analyse_Pairs <- PerSIMPER_onMatrix(Mat_Sampled, Group_Sampled, NomCluster = LETTERS[Unik_group], NS = FALSE,
      #                                    overall = FALSE)
      Analyse_Overall <- DNCI.ses_overall(Mat_Sampled, Group_Sampled, id = id, Nperm = Nperm, count = count)

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
