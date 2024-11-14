#' Per SIMPER function : identification of the main assembly process
#'
#'
#' This function is the basis of DNCImper package and DNCI analysis. It will permute the empirical matrix and produce the empirical as well as the randomized SIMPER profiles. Identify the main assembly process by comparing profiles. Permutations are fixed by rows, columns or both corresponding respectively, to niche, dispersal and niche+dispersal hypothesis.
#' The E index plot is produced to highlight the main assembly process. See Gibert & Escarguel 2019 Global Ecology and Biogeography for theory and more information on process identification.
#' More information in code and comments inside function file.
#' @param matrixSIMP Sample/Taxa matrix with sample in row and taxa in column
#' @param Groups Grouping vector, ex : c(1,1,1,1,2,2,2,2,2) : 2 groups only !!
#' @param count Display the number of permutation done, can be usefull with very large or small matrix, default = TRUE
#' @param dataTYPE Need to be set for presence/absence or abundance data ("count"), default = "prab" (presence_absence)
#' @param Nperm Number of permutation, default = 1000, should be change to 100 for robustness analysis
#' @param plotSIMPER Display the SIMPER, PerSIMPER and E index plots, default = TRUE
#' @param parallelComputing Run PerSIMPER on half of the available cores/nodes
#' @examples A <- DNCImper:::PerSIMPER(Matrix, Group)
#' @examples #where Matrix is a presence/absence matrix with taxa in column and sample in row
#' @examples #and Group is a vector with length() == number of rows/samples in Matrix, 2 groups ONLY
#' @examples #
#' @examples B <- DNCImper:::PerSIMPER(Matrix, Group, Nperm = 100, count = FALSE, plotSIMPER = FALSE)
#' @examples #In this example, same data are analysed, with 100 permutations, with no countdown and no plots
#' @importFrom graphics legend lines title
#' @importFrom stats median quantile sd
#' @importFrom utils combn
#' @importFrom vegan permatfull simper
#' @importFrom foreach %dopar% foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster detectCores stopCluster
#' @importFrom ggplot2 ggplot aes geom_boxplot scale_fill_manual ggtitle theme element_text scale_y_continuous labs scale_x_discrete
#'
#'


#
#############################################################
# Original Persimper function
# require packages: vegan, ggplot2
# by Corentin Gibert Bret
# original idea by Gilles Escarguel (LEHNA lab Lyon University)
# revised by Corentin Gibert Bret for parallel computing and infinite loop of permutation 2024-11-12
# revised by Jianjun Wang, 2019-02-01
# revised by Aurelien, 2019-07-19
#############################################################
#

PerSIMPER <- function(matrixSIMP,
                      Groups,
                      count = TRUE,   # default is true
                      dataTYPE = "prab",
                      Nperm=1000,
                      plotSIMPER = TRUE,
                      parallelComputing = FALSE){    # add the possibility to change the number of permutations

  #Package loading
  library(vegan)
  library(ggplot2)
  #library(dplyr)

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

  #matrixSIMP <- Stores the matrix to use in SIMPER analysis
  #(i.e. the presence/absence or the abundance distribution of taxa in at least 2 clusters of assemblages)
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

  ### Warning message ================
  if(parallelComputing == TRUE){
    library(doParallel)
    library(parallel)
    library(foreach)
    warning("Parallel computing is ON")
    warning("Count argument can't be used ")
    logMessage <- function(parallel_message){
      cat(parallel_message,"\n", file = "parallel_log.txt", append = TRUE)
    }


  }

  if(length(matrixSIMP[,1]) != length(Groups)){
    stop("Group vector and occurrence/abundance matrix are of different length")
  }

  if(Nperm >= 1000){
    warning("1000 or more permutation have been selected for the computation of each null model")
    warning("If you are using multigroup analysis and computation is slow, reduce Nperm and turn off plotSIMPER, turn on parallelComputing")
    if(parallelComputing == FALSE){
      warning("OPTION : Turn on parallel computing")}
  }

  warning("Per-SIMPER computation can be slow, especially if your presence/absence/abundance matrix is nearly empty or very large")

  ####### Start here =========================

  ###### Change abundance data (dataTYPE = "count") to relative abundance within sites with 1 as min(abundance) for the rarest taxa
  ###### This change is necessary A. for increasing the speed of swap algorithm, B. to increase comparability between datasets, C. because swap is impossible with non integer (< 1) values

  ### Count/abundance data :  ================

  if (dataTYPE == "count") {


    warning("WARNING : PER-SIMPER has been tested and mainly used on presence/absence data")
    warning("WARNING : You are using count/abundance data, please compare your results with pre/abs data")
    warning("WARNING : Overly abundant taxa can force PER-SIMPER toward very negative results")
    warning("WARNING : Relative abundance can lower this effect")

    if (min(matrixSIMP[1, matrixSIMP[1, ] != 0]) != 1) {
      print("Absolute abundance modified for relative abundance")
      tempMatrix <- matrixSIMP
      for (i in 1:length(matrixSIMP[, 1])) {
        diff0 <- which(matrixSIMP[i, ] != 0)
        min0 <- min(matrixSIMP[i, diff0])
        newLine <- round(matrixSIMP[i, diff0] * (1/min0))
        tempMatrix[i, diff0] <- newLine
      }
      matrixSIMP <- tempMatrix
    }
  }

  ### Presence/absence/occurrence data : =====

  ### SIMPER Analysis performed with vegan package
  AnaSimp <- vegan::simper(matrixSIMP, Groups)

  Contribution <- sort(AnaSimp[[1]]$average, decreasing = TRUE)
  #Replication in a vector (named 'Contribution') of the sorting
  #of species by their contribution to overall dissimilarity (OAD)

  Pourcent_Contribution <- ((Contribution)/sum(Contribution)) * 100
  #Conversion as a percentage of each species' contribution to the OAD

  ### Pritting a simper plot if arg = TRUE
  if (plotSIMPER == TRUE) {
    maxPlot <- max(Pourcent_Contribution) + 5 # Setting Y axis
    if (maxPlot > 100) {
      maxPlot <- 100
    }
    plot(Pourcent_Contribution, col = "brown2", type = "p",
         lwd = 1, ylab = "% contribution to dissimilarity",
         xlab = "Species", ylim = c(0, maxPlot))
  }

  ### Beginning of the permutation of the matrix ===============

  ### Permutation with fixed rows and columns sums
  dp2 <- vegan::permatfull(matrixSIMP, fixedmar = "both", mtype = dataTYPE,
                    times = Nperm)
  ### Permutation with fixed rows sums
  dp3 <- vegan::permatfull(matrixSIMP, fixedmar = "rows", mtype = dataTYPE,
                    times = Nperm)
  #Randomization of the matrixSIMP matrix ; permatfull need to be used in order to swap cells under
  #various conditions. permatswap only allow permutation with both fixed rows and fixed columns count.
  #mtype = "prab" is used for presence/absence data, this setting must be changed with "count" if abundance
  #analysis are performed.

  df2 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  df3 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  df4 <- matrix(nrow = Nperm, ncol = length(matrixSIMP[2,]))
  #Generating matrices that will store the results (the ranked contribution of species to the OAD)
  #of the 1000 permutations of the original matrix

  ### Permutation for fixed columns sums can be problematic (stuck in infinite loop,
  ### looking for new permutations)
  ### This permutation is decomposed to check for infinite loop
  ### The algorithm can switch to manual swapping to solve infinite looping

  if(parallelComputing == TRUE){
    warning("If permutation end up in infinite loop, log will be written in WD")

    #Initiation of CPU nodes
    cl <- parallel::makeCluster(parallel::detectCores() / 2) # use one less core than the available to not overload the system registerDoParallel(cl
    doParallel::registerDoParallel(cl)
    ParallelPerSIMPER <- foreach::foreach(i = 1:Nperm) %dopar% {

      SWAPcount <- 0

      ### Permutation with fixed columns will repeat until no rows or columns are empty
      repeat {

        ### The algorithm will try 200 times to permute before manual override
        SWAPcount <- SWAPcount + 1
        ### if v == T, repeat will break
        v <- T
        dp4 <- vegan::permatfull(matrixSIMP, fixedmar = "columns",
                                  mtype = dataTYPE, times = 1)
        ### Check if rows are empty (can happen with this algorithm)
        ### If one or more is empty, v <- F and the repeat loop will continue
        if (length(which(apply(dp4$perm[[1]], 1, sum) ==
                         0)) != 0) {
          v <- F
        }

        ### If this repeat fail 200 times, manual override is done
        if (v == F & SWAPcount > 200) {

          logMessage(paste(warning("Permutated rows are still empty: manual override"), "\n"))

          repeat {
            for (j in 1:length(dp4$perm[[1]][, 2])) {
              ### Finding the empty rows
              if (sum(dp4$perm[[1]][j, ]) == 0) {
                v <- FALSE
                ### Columns counts and selection of the richest quartile of columns
                ### Then random selection of a cells in this column

                tempColSum <- apply(dp4$perm[[1]], 2,
                                    sum)
                tempHigh_Col <- sample(which(apply(dp4$perm[[1]],
                                                   2, sum) >= quantile(tempColSum)[4]), 1)
                tempCel <- which(dp4$perm[[1]][, tempHigh_Col] > 0)
                SelectedRows <- tempCel[which(apply(dp4$perm[[1]][tempCel,
                ], 1, sum) >= quantile(apply(dp4$perm[[1]][tempCel,
                ], 1, sum))[4])]

                SelectedCell <- sample(SelectedRows, 1)
                dp4$perm[[1]][SelectedCell, tempHigh_Col] <- 0

                ### Within one of the richest column, a 1 is moved to the empty row j selected previously in the loop
                dp4$perm[[1]][j, tempHigh_Col] <- 1

                ### Check if all rows are empty, if not V = TRUE and repeat are stopped by break instructions
                if (length(which(apply(dp4$perm[[1]],
                                       1, sum) == 0)) == 0) {
                  v <- TRUE
                }
              }
            }
            if (v == TRUE) {
              break
            }
          }
        }

        ### Check if all rows are empty in the first repeat loop (not the manual override)
        if (length(which(apply(dp4$perm[[1]], 1, sum) ==
                         0)) == 0) {
          break
        }
      }

      ### SIMPER analysis are repeated on permuted matrix via vegan package
      ### Results are sorted by contribution as in the SIMPER analysis of the empirical dataset (AnaSimp)
      simp2 <- vegan::simper(dp2$perm[[i]], Groups)
      simp3 <- vegan::simper(dp3$perm[[i]], Groups)
      simp4 <- vegan::simper(dp4$perm[[1]], Groups)

      list(simp2 = sort(simp2[[1]]$average, decreasing = TRUE),
           simp3 = sort(simp3[[1]]$average, decreasing = TRUE),
           simp4 = sort(simp4[[1]]$average, decreasing = TRUE))

    }

    ### doparallel output are list, the sorting of permutation results are adapted to this new format
    for(i in 1:Nperm)
    {
      df2[i,] <- (ParallelPerSIMPER[[i]]$simp2/sum(ParallelPerSIMPER[[i]]$simp2)) * 100
      df3[i,] <- (ParallelPerSIMPER[[i]]$simp3/sum(ParallelPerSIMPER[[i]]$simp3)) * 100
      df4[i,] <- (ParallelPerSIMPER[[i]]$simp4/sum(ParallelPerSIMPER[[i]]$simp4)) * 100
    }

    ### Stopping parallel computing
    stopCluster(cl)
  }

  ### Start of the sequential computing, this section is ignored if parallelComputing == TRUE
  if(parallelComputing == FALSE){

    for (i in 1:Nperm) {

      ### Simple Count to track Per-SIMPER speed
      if (count == TRUE && i < 100 || count == TRUE && i >
          round(Nperm * 0.9)) {
        print(i)
      }
      SWAPcount <- 0

      ### Permutation with fixed columns will repeat until no rows or columns are empty
      repeat {

        ### The algorithm will try 200 times to permut before manual override
        SWAPcount <- SWAPcount + 1
        ### if v == T, repeat will break
        v <- T
        dp4 <- permatfull(matrixSIMP, fixedmar = "columns",
                          mtype = dataTYPE, times = 1)
        ### Check if rows are empty (can happen with this algorithm)
        ### If one or more is empty, v <- F and the repeat loop will continue
        if (length(which(apply(dp4$perm[[1]], 1, sum) ==
                         0)) != 0) {
          v <- F
        }

        ### If this repeat fail 200 times, manual override is done
        if (v == F & SWAPcount > 200) {
          repeat {
            for (j in 1:length(dp4$perm[[1]][, 2])) {
              ### Finding the empty rows
              if (sum(dp4$perm[[1]][j, ]) == 0) {
                v <- FALSE
                ### Columns counts and selection of the richest quartile of columns
                ### Then random selection of a cells in this column

                tempColSum <- apply(dp4$perm[[1]], 2,
                                    sum)
                tempHigh_Col <- sample(which(apply(dp4$perm[[1]],
                                                   2, sum) >= quantile(tempColSum)[4]), 1)
                tempCel <- which(dp4$perm[[1]][, tempHigh_Col] > 0)
                SelectedRows <- tempCel[which(apply(dp4$perm[[1]][tempCel,
                ], 1, sum) >= quantile(apply(dp4$perm[[1]][tempCel,
                ], 1, sum))[4])]

                SelectedCell <- sample(SelectedRows, 1)
                dp4$perm[[1]][SelectedCell, tempHigh_Col] <- 0

                ### Within one of the richest column, a 1 is moved to the empty row j selected previously in the loop
                dp4$perm[[1]][j, tempHigh_Col] <- 1

                ### Check if all rows are empty, if not V = TRUE and repeat are stopped by break instructions
                if (length(which(apply(dp4$perm[[1]],
                                       1, sum) == 0)) == 0) {
                  v <- TRUE
                }
              }
            }
            if (v == TRUE) {
              break
            }
          }
        }

        ### Check if all rows are empty in the first repeat loop (not the manual override)
        if (length(which(apply(dp4$perm[[1]], 1, sum) ==
                         0)) == 0) {
          break
        }
      }

      ### SIMPER analysis are repeated on permuted matrix via vegan package
      ### Results are sorted by contribution as in the SIMPER analysis of the empirical dataset (AnaSimp)
      simp2 <- simper(dp2$perm[[i]], Groups)
      simp3 <- simper(dp3$perm[[i]], Groups)
      simp4 <- simper(dp4$perm[[1]], Groups)
      df2[i, ] <- sort(simp2[[1]]$average, decreasing = TRUE)
      df3[i, ] <- sort(simp3[[1]]$average, decreasing = TRUE)
      df4[i, ] <- sort(simp4[[1]]$average, decreasing = TRUE)
      df2[i, ] <- (df2[i, ]/sum(df2[i, ])) * 100
      df3[i, ] <- (df3[i, ]/sum(df3[i, ])) * 100
      df4[i, ] <- (df4[i, ]/sum(df4[i, ])) * 100
    }
  }

  ### The Nperm SIMPER results are sorted in the result matrix by rows in order to
  ### produce envelope for PER-SIMPER; e.g. the Nperm/2 ==> median results
  dn2 <- apply(df2, 2, sort)
  dn3 <- apply(df3, 2, sort)
  dn4 <- apply(df4, 2, sort)

  ### Addition of the three null models envelope to the SIMPER empirical plot created Line 60
  if (plotSIMPER == TRUE) {
    lines(dn3[0.975 * Nperm, ], lty = "dotted", lwd = 2,
          col = "#009900")
    lines(dn3[0.025 * Nperm, ], lty = "dotted", lwd = 2,
          col = "#009900")
    lines(dn4[0.975 * Nperm, ], lty = "dotted", lwd = 2,
          col = "#FF9933")
    lines(dn4[0.025 * Nperm, ], lty = "dotted", lwd = 2,
          col = "#FF9933")
    title("SIMPER (in red) and PER-SIMPER profiles")
    legend(x = "topright", bty = "n", legend = c("SIMPER profil",
                                                 "Rows fixed", "Col fixed"), col = c("brown2", "#009900",
                                                                                     "#FF9933"), pch = c(15, 15, 15))
  }


  ###########################################################################
  #### The following is the calculation and the illustration of E index ####
  ####  E = Log of the sum of square deviations with empirical profile  ####
  ###########################################################################

  up <- 0.975 * Nperm
  lo <- 0.025 * Nperm
  med <- 0.5 * Nperm
  obs <- Pourcent_Contribution

  ### Orange is column-fixed null model (i.e. dispersal)
  ### Green is row-fixed null model (i.e. niche)
  ### Blue is rows+columns-fixed null model (i.e. niche and dispersal)
  Orange <- dn4
  Blue <- dn2
  Green <- dn3
  # Ranked % of contribution to OAD of empirical and simulated profiles

  ### Prepping vector for storing sum of log squared deviation from SIMPER empirical profile
  VectorEcartCarreOrangeLog <- vector(mode = "numeric", Nperm)
  VectorEcartCarreGreenLog <- vector(mode = "numeric", Nperm)
  VectorEcartCarreBlueLog <- vector(mode = "numeric", Nperm)
  for (i in 1:Nperm) {

    ### Prepping vector for storing simple deviation from SIMPER empirical profile with the 3 nulls profiles
    SommeEcartCarreOrange <- vector(mode = "numeric", length = length(Orange[1,
    ]))
    SommeEcartCarreGreen <- vector(mode = "numeric", length = length(Green[1,
    ]))
    SommeEcartCarreBlue <- vector(mode = "numeric", length = length(Blue[1,
    ]))

    ### Computation of the squared deviation from empirical profile
    ### Squared deviation ensure that all deviation are positives
    for (j in 1:length(obs)) {
      SommeEcartCarreOrange[j] <- (Orange[i, j] - obs[j])^2
      SommeEcartCarreGreen[j] <- (Green[i, j] - obs[j])^2
      SommeEcartCarreBlue[j] <- (Blue[i, j] - obs[j])^2
    }

    ### Sum of the deviation are transformed into log
    VectorEcartCarreOrangeLog[i] <- log10(sum(SommeEcartCarreOrange))
    VectorEcartCarreGreenLog[i] <- log10(sum(SommeEcartCarreGreen))
    VectorEcartCarreBlueLog[i] <- log10(sum(SommeEcartCarreBlue))
  }

  ### The mean of all these deviation are stored
  meanCarreOrangeLog <- mean(VectorEcartCarreOrangeLog)
  meanCarreGreenLog <- mean(VectorEcartCarreGreenLog)
  meanCarreBlueLog <- mean(VectorEcartCarreBlueLog)
  DataMeanCarreLog <- data.frame(Orange = VectorEcartCarreOrangeLog,
                                 Blue = VectorEcartCarreBlueLog, Green = VectorEcartCarreGreenLog)


  ### Plotting of permuted null models
  ## BOXPLOT with 95 % intervals of E (Log of the sum of square deviations with empirical profile) index
  if (plotSIMPER == TRUE) {
    Ax <- c("Fixed columns", "Both fixed", "Fixed rows")
    y <- DataMeanCarreLog
    df <- data.frame(Permutation_model = Ax, y0 = quantile(y$Orange,
                                                           0.025), y25 = quantile(y$Orange, 0.25), y50 = median(y$Orange),
                     y75 = quantile(y$Orange, 0.75), y100 = quantile(y$Orange, 0.975))
    df[2, 2] = quantile(y$Blue, 0.025)
    df[2, 3] = quantile(y$Blue, 0.25)
    df[2, 4] = median(y$Blue)
    df[2, 5] = quantile(y$Blue, 0.75)
    df[2, 6] = quantile(y$Blue, 0.975)
    df[3, 2] = quantile(y$Green, 0.025)
    df[3, 3] = quantile(y$Green, 0.25)
    df[3, 4] = median(y$Green)
    df[3, 5] = quantile(y$Green, 0.75)
    df[3, 6] = quantile(y$Green, 0.975)
    #print(df)
    E <- ggplot(df, aes(Permutation_model, fill = Permutation_model)) +
      geom_boxplot(aes(ymin = y0, lower = y25, middle = y50,
                       upper = y75, ymax = y100), stat = "identity") +
      scale_fill_manual(values = c("#00CCFF", "#FF9933",
                                   "#009900")) + ggtitle("The lower E, the closer the simulated profile to empirical SIMPER profile") +
      theme(plot.title = element_text(lineheight = 2)) +
      scale_y_continuous(name = " E (Log of the sum of square deviations with empirical profile)") +
      labs(fill = "Permutation models") + scale_x_discrete(name = "Permutation models")
    print(E)
  }

  ### Storing all results and PerSIMPER plots parameters in a list() =============
  ListResults <- list(EcartCarreLog = DataMeanCarreLog, mat = matrixSIMP,
                      ContriPercentage = Pourcent_Contribution, UpOrange = dn4[up,
                      ], DownOrange = dn4[lo, ], MedOrange = dn4[med,
                      ], UpBlue = dn2[up, ], DownBlue = dn2[lo, ], MedBlue = dn2[med,
                      ], UpGreen = dn3[lo, ], DownGreen = dn3[up, ], MedGreen = dn3[med,
                      ], dnOrange = dn4, dnBlue = dn2, dnGreen = dn3)

  return(ListResults)

}
