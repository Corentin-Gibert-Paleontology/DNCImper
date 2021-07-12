#' Per SIMPER function : identification of the main assembly process
#'
#'
#' This function is the basis of DNCImper package and DNCI analysis. It will permute the empirical matrix and produce the empirical as well as the randomized SIMPER profiles. Identify the main assembly process by comparing profiles. Permutations are fixed by rows, columns or both corresponding respectively, to niche, dispersal and niche+dispersal hypothesis.
#' The E index plot is produced to highlight the main assembly process. See Gibert & Escarguel 2019 Global Ecology and Biogeography for theory and more information on process identification.
#' More information in code and comments inside function file.
#' @param matrixSIMP Sample/Taxa matrix with sample in row and taxa in column
#' @param Groups Grouping vector, ex : c(1,1,1,1,2,2,2,2,2) : 2 groups only !!
#' @param leg Display the legend on PerSIMPER profile, default = FALSE
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataType Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @examples A <- DNCImper:::PerSIMPER(Matrix, Group)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 2 groups ONLY
#' @examples #
#' @examples B <- DNCImper:::PerSIMPER(Matrix, Group, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
#'
#'


#
#############################################################
# Original Persimper function
# require packages: vegan, ggplot2
# by Corentin Gibert
# original idea by Gilles Escarguel (LEHNA lab Lyon University)
# revised by Jianjun Wang, 2019-02-01
# revised by Aurelien, 2019-07-19
#############################################################
#

PerSIMPER <- function(matrixSIMP,
                          Groups,
                          leg = FALSE,
                          count = TRUE,   # default is true
                          dataTYPE = "prab",
                          Nperm=1000,
                          plotSIMPER = TRUE){    # add the possibility to change the number of permutations

library(vegan)
library(ggplot2)
library(dplyr)
  ################################################################################
  ## for every problems or questions, please contact me by mail or ResearchGate: ##
  ##### at corentingibert@gmail.com | corentin.gibert@univ-poitiers.fr ############
  ########### https://www.researchgate.net/profile/Corentin_Gibert ################
  #################################################################################

  #### GitHub : https://github.com/Corentin-Gibert-Paleontology

  # The PerSIMPER function requests exactly the same required arguments (matrixSIMP & Groups)
  # as the functions used to compute SIMPER method in R, i.e. a presence/absence matrix and
  # a vector encoding cluster information.

  #Arguments :

  #matrixSIMP <- Stores the matrix to use in SIMPER analysis (i.e. the result
  #of the presence/absence or the abundance distribution of taxa in at least 2 clusters of assemblies)
  # LOCALITIES in LINES
  # TAXA in COLUMNS

  #Groups <- A vector allowing to assign to the lines of the SIMPER matrix (matrixSIMP) a cluster of
  #localities ; for example in a matrix of 10 lines built from 2 sets of localities
  #of the same size : 1, 1, 1, 1, 1, 2, 2, 2, 2, 2. Strings are accepted, e.g., RegionA, Region, A, Region B, etc...

  #log <- the "log" argument allows to log the y-axis (i.e. the percentage contribution to the OAD of the species)

  #leg <- the "leg" argument allows a legend to be displayed on the Per-SIMPER profile

  #count <- the "count" argument allows a Screen output of the number of iterations performed.
  #This option is used to indicate if the permutation function is unable to swap the matrix cells.
  #This incapacity is usually the result of a matrix too sparse in data (too many cells at 0).

  #dataTYPE <- the "dataTYPE" argument allows to choose between presence/absence data or abundance data. By default
  #the algorithm use presence/absence permutation for presence/absence data, if you use abundance dataset, you need
  #to write "count" in the dataTYPE argument.
  #e.g. dataTYPE = "count"

  #Nperm <- number of matrix permutation

  AnaSimp <- simper(matrixSIMP, Groups)   # summary(AnaSimp)
  #Classical SIMPER analysis computed on the compared groups

  Contribution <- sort(AnaSimp[[1]]$average, decreasing = TRUE)
  #Replication in a vector (named 'Contribution') of the sorting
  #of species by their contribution to overall dissimilarity (OAD)

  Pourcent_Contribution <- ((Contribution)/sum(Contribution))*100
  #Conversion as a percentage of each species' contribution to the OAD

  if(plotSIMPER == TRUE)
  {
   # if(logSIMPER == TRUE) {
   # plot(Pourcent_Contribution, col = "brown2", log="y", type="p",lwd = 1.5,
   #      ylab ="%  contribution to dissimilarity", xlab="Species")
   # } else {
     plot(Pourcent_Contribution, col = "brown2", type="p",lwd = 1,
         ylab ="% contribution to dissimilarity", xlab="Species")
  #Ploting SIMPER results in percentage (in Log or not)
  }


  # set.seed(123456)  # add by jjwang
  dp2 <- permatfull(matrixSIMP, fixedmar = "both", mtype = dataTYPE,  times = Nperm) #prab
  # set.seed(123456)  # add by jjwang
  dp3 <- permatfull(matrixSIMP, fixedmar ="rows" , mtype = dataTYPE, times = Nperm)  #prab

  #Randomization of the matrixSIMP matrix ; permatfull need to be used in order to swap cells under
  #various conditions. permatswap only allow permutation with both fixed rows and fixed columns count.
  #mtype = "prab" is used for presence/absence data, this setting must be changed with "count" if abundance
  #analysis are performed.

  df2 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  df3 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  df4 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  #Generating matrices that will store the results (the ranked contribution of species to the OAD)
  #of the 1000 permutations of the original matrix

  jj = 0

  for (i in 1:Nperm)  {
    if(count == TRUE && i < 100 || count == TRUE && i > round(Nperm*0.90) ){print(i)}
    #Screen output of the number of iterations performed.
    #This option is used to indicate if the permutation function is unable to swap the matrix cells.
    #This incapacity is usually the result of a matrix too sparse in data (too many cells at 0).
    repeat {
      v <- T
      dp4 <- permatfull(matrixSIMP, fixedmar = "columns", mtype = dataTYPE, times = 1)  #prab
      for(j in 1:length(dp4$perm[[1]][,2])) {
        if(sum(dp4$perm[[1]][j,]) == 0){v <- FALSE
             
        if(SWAPcount > 50)
          {
          ### Looking for cells to swap from rich taxa and rich locality to empty locality
          tempColSum <- apply(dp4$perm[[1]], 2, sum)
          tempHigh_Col <- which(apply(dp4$perm[[1]], 2, sum) > median(tempColSum))
          tempMoove <- sample(tempHigh_Col, 1)
          tempColSel <- dp4$perm[[1]][,tempMoove]
          tempCel <-  which(tempColSel  == 1)
          tempHigh_Row <- which(apply(dp4$perm[[1]][tempCel,], 1, sum) > median(apply(dp4$perm[[1]][tempCel,], 1, sum))) 
          tempMoove2 <- sample(tempHigh_Row, 1)
          
          #### Swapping
          dp4$perm[[1]][tempCel[tempMoove2], tempMoove] <- 0
          dp4$perm[[1]][j,tempMoove] <- 1
          v <- TRUE
          }                                  
                                       
         }
        }
      if(v == TRUE) break
      }

    simp2 <- simper(dp2$perm[[i]], Groups)
    simp3 <- simper(dp3$perm[[i]], Groups)
    simp4 <- simper(dp4$perm[[1]], Groups)
    #SIMPER analysis performed on each permutated matrix

    df2[i,] <- sort(simp2[[1]]$average, decreasing = TRUE)
    df3[i,] <- sort(simp3[[1]]$average, decreasing = TRUE)
    df4[i,] <- sort(simp4[[1]]$average, decreasing = TRUE)
    #Storage of SIMPER results (ranked contribution to OAD)

    df2[i,] <- (df2[i,]/sum(df2[i,]))*100
    df3[i,] <- (df3[i,]/sum(df3[i,]))*100
    df4[i,] <- (df4[i,]/sum(df4[i,]))*100
    #Conversion to percentage of SIMPER results
  }   #

  dn2 <- apply(df2, 2, sort)
  dn3 <- apply(df3, 2, sort)
  dn4 <- apply(df4, 2, sort)

if(plotSIMPER == TRUE){
lines(dn3[0.975*Nperm,], lty="dotted", lwd=2, col="cornflowerblue")
lines(dn3[0.025*Nperm,], lty="dotted", lwd=2, col="cornflowerblue")
lines(dn4[0.975*Nperm,], lty="dotted", lwd=2, col="orange2")
lines(dn4[0.025*Nperm,], lty="dotted", lwd=2, col="orange2")
#Plot of the upper and lower limit of the confidence intervals for fixed rows,
#fixed columns and fixed rows and columns
title("SIMPER (in red) and PER-SIMPER profiles")

legend(x="topright", bty="n",legend=c("SIMPER profil", "Rows fixed", "Col fixed"),
  col=c("brown2", "cornflowerblue","orange2"), pch=c(15,15,15))
  }

  up <- 0.975*Nperm # ex for 100 permutations it will used 97
  lo <- 0.025*Nperm # ex for 100 permutations it will used 2
  med <- 0.5*Nperm

  ###########################################################################
  #### The following is the calculation and the illustration of E index ####
  ####  E = Log of the sum of square deviations with empirical profile  ####
  ###########################################################################

  obs <- Pourcent_Contribution
  Orange <- dn4
  Blue <- dn2
  Green <- dn3
  # Ranked % of contribution to OAD of empirical and simulated profiles

  VectorEcartCarreOrangeLog <- vector(mode = "numeric", Nperm)
  VectorEcartCarreGreenLog <- vector(mode = "numeric", Nperm)
  VectorEcartCarreBlueLog <- vector(mode = "numeric", Nperm)

  for(i in 1:Nperm)  {
    SommeEcartCarreOrange <- vector(mode = "numeric", length = length(Orange[1,]))
    SommeEcartCarreGreen <- vector(mode = "numeric", length = length(Green[1,]))
    SommeEcartCarreBlue <- vector(mode = "numeric", length = length(Blue[1,]))

    for(j in 1:length(obs))    {
      SommeEcartCarreOrange[j] <-  (Orange[i,j] - obs[j])^2
      SommeEcartCarreGreen[j] <-   (Green[i,j] - obs[j])^2
      SommeEcartCarreBlue[j] <-    (Blue[i,j] - obs[j])^2
    }
    # Computation of square deviations with empirical profile (obs)

    #Mise en log des carre des ecarts pour symetriser la distribution
    VectorEcartCarreOrangeLog[i] <- log10(sum(SommeEcartCarreOrange))
    VectorEcartCarreGreenLog[i] <- log10(sum(SommeEcartCarreGreen))
    VectorEcartCarreBlueLog[i] <- log10(sum(SommeEcartCarreBlue))
  }
  # Log conversion of the sum of square deviations

  meanCarreOrangeLog <- mean(VectorEcartCarreOrangeLog)
  meanCarreGreenLog <- mean(VectorEcartCarreGreenLog)
  meanCarreBlueLog <- mean(VectorEcartCarreBlueLog)

  DataMeanCarreLog <- data.frame(Orange=VectorEcartCarreOrangeLog,
                                 Blue=VectorEcartCarreBlueLog,
                                 Green=VectorEcartCarreGreenLog)

if(plotSIMPER == TRUE)
{

   ## BOXPLOT with 95 % intervals of E (Log of the sum of square deviations with empirical profile) index

   Ax <- c("Fixed columns", "Both fixed", "Fixed rows")
   y <- DataMeanCarreLog
   df <- data.frame(
     Permutation_model = Ax,
     y0 = quantile(y$Orange, 0.025),
     y25 = quantile(y$Orange, 0.25),
     y50 = median(y$Orange),
     y75 = quantile(y$Orange, 0.75),
     y100 = quantile(y$Orange, 0.975)
   )
   df[2,2] = quantile(y$Blue, 0.025)
   df[2,3] = quantile(y$Blue, 0.25)
   df[2,4] = median(y$Blue)
   df[2,5] = quantile(y$Blue, 0.75)
   df[2,6] = quantile(y$Blue, 0.975)

   df[3,2] = quantile(y$Green, 0.025)
   df[3,3] = quantile(y$Green, 0.25)
   df[3,4] = median(y$Green)
   df[3,5] = quantile(y$Green, 0.75)
   df[3,6] = quantile(y$Green, 0.975)
   # Extraction of quantiles of interest

  # comment out by jjwang
   print(df)   # comment out by jjwang, 2019-01-30

    E <- ggplot(df, aes(Permutation_model, fill = Permutation_model)) +
      geom_boxplot(
        aes(ymin = y0, lower = y25, middle = y50, upper = y75, ymax = y100),
        stat = "identity") + scale_fill_manual(values=c("#CCCCCC", "#FF6600", "#00CCFF")) +
      ggtitle("The lower E, the closer the simulated profile to empirical SIMPER profile") +
      theme(plot.title = element_text(lineheight= 2)) +
      scale_y_continuous(name=" E (Log of the sum of square deviations with empirical profile)")+
      labs(fill="Permutation models")+
      scale_x_discrete(name = "Permutation models")

    print(E)

  #Computation of ggplots2 boxplot for E index

  # outputs
}

  ListResults <- list(EcartCarreLog = DataMeanCarreLog,
                      mat = matrixSIMP, ContriPercentage = Pourcent_Contribution,
                      UpOrange = dn4[up,], DownOrange = dn4[lo,], MedOrange = dn4[med,],
                      UpBlue = dn2[up,], DownBlue = dn2[lo,], MedBlue = dn2[med,],
                      UpGreen = dn3[lo,], DownGreen = dn3[up,], MedGreen = dn3[med,],
                      dnOrange = dn4, dnBlue = dn2, dnGreen = dn3)

  return(ListResults)
}



