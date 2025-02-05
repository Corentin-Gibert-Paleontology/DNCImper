# PER-SIMPER-DNCI functions package: V 1.7

A bug has been removed from core function PerSIMPER and PerSIMPER_overall that could have affected a bit the DNCI results. 

#### Warning -- Use R >= 4.0 <br/>
If the error message "can't find function X" appears:<br/>
use DNCImper:::NAME_Function<br/> 
e.g. DNCImper:::PerSIMPER(Matrix, Group)<br/> 

To access informations on functions, you can use the classical ? call in R/RStudio<br/> 
e.g. ?DNCImper:::PerSIMPER())<br/>
A large documentation, argument explanation and examples are associated with each functions<br/>

#######################################<br/>
HI and THANK YOU ########################<br/>
FOR YOUR INTEREST IN PER-SIMPER ###########<br/>
 AND DISPERSAL-NICHE-CONTINUUM INDEX #####<br/>
#######################################<br/>
BY CORENTIN GIBERT (author, creator) ##########<br/> 
AND GILLES ESCARGUEL (author) ##############<br/>
ANNIKA VILMI AND JIANJUN WANG ###########<br/>
Thanks to MAXIME LOPEZ for bug correction #####<br/>
#######################################<br/>

WELCOME to the DNCImper PACKAGE README.txt<br/>

The main function to compute DNCI is DNCI_multigroup()<br/>
It requires 2 packages (vegan, ggplot2)<br/>

### Install by using (devtools package) and write : devtools::install_github("Corentin-Gibert-Paleontology/DNCImper")<br/>
(error messages like "ERROR 404 not found" or "can't get to the repo" or "this R is version 3.XXXX, package 'DNCImper' requires R >= 3.10" should be solved with updating your R. to the new 4.0 version)<br/>

### If devtools can't download the package, download manually the GitHub repository and install it with devtools:::install(DNCImper). Set your Working Directory in the folder containing the unzipped DNCImper repository beforehand.<br/>
Second option is to manually install the DNCImper[...].tar.gz file from the lastest release (can be downloaded on the right side of my GitHub account) via the package window of RStudio. 

You can find a tutorial (from v1.5) named "Tutorial_for_R.R" in https://github.com/Corentin-Gibert-Paleontology/PER-SIMPER-DNCI_Tutorial<br/>

This file will help you to use PerSIMPER,<br/> 
DNCI and other functions associated with<br/>
Vilmi, Gibert et al. 2021 and Gibert & Escarguel<br/> 
2019. 

Please feel free to contact me to highlight<br/>
bugs, add suggestions or ask for helps.<br/>

Contact:<br/> 
corentingibert@gmail.com<br/> 
or<br/>
corentin.gibert-bret@univ-lille.fr<br/> 
annika.vilmi@gmail.com<br/>
annika.vilmi@syke.fi<br/>

ResearchGate:<br/>
https://www.researchgate.net/profile/Corentin_Gibert<br/>
https://www.researchgate.net/profile/Annika_Vilmi<br/>

GitHub:<br/> 
https://github.com/Corentin-Gibert-Paleontology


