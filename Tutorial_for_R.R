##############################################################
###### EXAMPLE FILE for PerSIMPER and DNCI function usage ####
##############################################################

### Loading library, functions, datas & working directory
library(vegan)
library(plyr)
library(ggplot2)

setwd("Your working directory")

load("Data_Matrix_Tutorial.RData")
load("Data_Group_Turorial.RData")

################################################################################
##### PerSIMPER usage for 2 groups #############################################
################################################################################
source("PerSIMPER.R") #Explanation for arguments inside function files
X1 <- PerSIMPER(Matrix, Group) #Default = 1000 permutations, Plot = TRUE

#Produce 2 plots and a list() with plots data. A simper plot + PerSIMPER profiles, and the E index plot.

################################################################################
##### DNCI - Dispersal-Niche Continuum Index usage for 2 groups ################
################################################################################
source("PerSIMPER.R")
source("DNCI_ses.R") #Explanation for arguments inside function files
source("DNCI_multigroup.R") #Explanation for arguments inside function files

X2 <- DNCI.ses(Matrix, Group) #Original DNCI function, limited to 2 groups
or
X2 <- DNCI_multigroup(Matrix, Group, plotSIMPER = TRUE) #Wraper, 2 or more groups
#Produce 2 plots -if plotSIMPER = TRUE- and a data.frame() with the DNCI result, 
#                                       associated with its confidence interval.

######
###### Results can be made more robust by making group even ## REPEAT X times ## 
X2 <- DNCI_multigroup(Matrix, Group, plotSIMPER = TRUE, symmetrize = TRUE)

################################################################################
###### DNCI usage for more than 2 groups #######################################
################################################################################
source("PerSIMPER.R")
source("DNCI_ses.R")
source("DNCI_multigroup.R")
load("Data_Matrix_Tutorial_4Groups.RData")
load("Data_Group4_Tutorial.RData")

X3 <- DNCI_multigroup(Matrix_4groups, Group4)
#Produce 2 plots by pairs of groups -if plotSIMPER = TRUE- and a data.frame() 
#                                           with the DNCI result for each pairs, 
#                                           associated with their respective confidence interval.
######
###### Results can be made more robust by making pairs of groups even ## REPEAT X times ## 
X3 <- DNCI_multigroup(Matrix_4groups, Group4, plotSIMPER = TRUE, symmetrize = TRUE)

######################## OVERALL DNCI - PerSIMPER ##############################
############## COMPUTATION ON THE OVERALL DISSIMILARITY ########################
##################### BETWEEN MORE THAN 2 GROUPS ###############################
############## THIS ANALYSIS IS NOT BY PAIRS OF GROUPS #########################
################################################################################
source("PerSIMPER_Overall.R") #Explanation for arguments inside function files
source("DNCI_ses_overall.R")  #Explanation for arguments inside function files
source("DNCI_ses_overall_symmetrized.R")
load("Data_Matrix_Tutorial_4Groups.RData")
load("Data_Group4_Tutorial.RData")

##### PerSIMPER Overall #### QUALITATIVE #######################################
X4 <- PerSIMPER_overall(Matrix_4groups, Group4)

##### DNCI Overall ##### QUANTITATIVE ##########################################
X5 <- DNCI.ses_overall(Matrix_4groups, Group4)
#Produce 2 plots -if plotSIMPER = TRUE- and a data.frame() with the overall DNCI result,
#                                       associated with its confidence interval.

######
###### Results can be made more robust by making group even ## REPEAT X times ## 
X6 <- DNCI.ses_overall_symmetrized(Matrix_4groups, Group4, NbrReRun = 10)
###### NbrReRun = argument to choose the number of iterations for resampling ###
