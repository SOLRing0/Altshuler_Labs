##-------- Main Script -----------
#This code will output all aerodynamic, lab manipulation and camber data
working_directory = "ENTERYOURDIRECTORYHERE"
#------- Run the Initialize (definitions of labels/themes/functions) -------
setwd(working_directory) 
source('InitializePlotting.R')
#------- Run all Aerodynamic data processing -------
source('AerodynamicData.R')
#------- Run the Camber-izer -------
source('Camber-izer.R')
#------- Run the Stats ------
source('Stats.R')
#------- Run the GeoMorph-izer -------
source('GeoMorph_projection-izer.R')
#Opening Order:  1 - Labman; 2- Invivo ; 3 - WT

##### Create Wireframe IF NEEDED ##### 
mean_shape <- mshape(labman.3D) #finds mean shape of the lab manipulations 
linked_landmarks <- define.links(mean_shape)


