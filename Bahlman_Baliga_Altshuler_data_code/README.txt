Data, code, and 3D printing files for:
	"Flight muscle power increases with strain amplitude and decreases 
	with cycle frequency for small birds" by JW Bahlman, VB Baliga, and 
	DL Altshuler

Questions/comments on Figshare contents: 
	Please contact Vikram Baliga (vbaliga@zoology.ubc.ca)

To re-create data analyses and figures, use the R scripts:
	1) Each .R file in the root directory is numbered. Please run them 
	   in order. Objects & functions in lower-numbered scripts may be 
	   carried into subsequent scripts. All raw data, which are located 
	   in the "raw_data" folder are imported and analyzed using the .R
	   scripts.
	   
	2) File "01_all_functions.R" has custom-written code for the import 
	   and analysis of muscle physiology data in R. You may simply 
	   `source()` this file and move on. All of this code was later re-
	   written as the R package 'workloopR', which was subsequently 
	   peer-reviewed by the rOpenSci Project and published in JOSS. 
	   
		  Code for workloopR is availabe at: 
		  https://github.com/ropensci/workloopR
		  
		  An accompanying paper is available at: 
		  https://doi.org/10.21105/joss.01856 
	   
	3) File "02_data_import.R" then imports and analyzes raw data from 
	   each bird. Running this script ultimately creates one mega-table 
	   of all relevant summary data from experiments (the object 
	   `allsummarydat` in the code).

	4) File "03_stats.R" runs our statistical analysis, which includes 
	   the competing of mixed-models and checks of model fit and 
	   sufficiency.

	5) File "04_figure_generator.R" provides the plots that ultimately 
	   make up all the figures in both the main text and in figure 
	   supplements.
	
Folder - Brace part piles:
	This folder contains the part files for duplicating the brace that 
	was used to immobilize the bird’s thorax during experiments. The 
	brace is composed of two 3D printed pieces (top and bottom), which 
	are connected with 4 bolts. The folder contains files with the 
	extension .STL which can be used to 3D print your own copies of 
	the brace. The folder also contains files with extension .FCstd, 
	which are for the software FreeCAD and can be used to edit the 
	brace design.
	
File - AllBirds_LoopUpTable.xlsx:
	An Excel file that helps match individual birds to experiment dates
	and specific experiments that the bird was used for. This file 
	itself is not be necessary to re-create analyses, but rather is 
	included to help you understand how the raw data are organized.
	