Data and code for:
	"Specializations in optic flow encoding in the pretectum of hummingbirds
	and zebra finches" by Smythe G, Baliga V, Gaede AH, Wylie DR, and Altshuler
	DL

Any questions about this repository can be directed to:
	Vikram B. Baliga (vbaliga@zoology.ubc.ca)

How this is organized:
	All analytic work and primary construction of figures were performed using
	R 4.1.2 in Rstudio 2021.09.1. An .Rproj file is included in the root 
	directory, and its use is strongly recommended to help set up working
	directories etc in Rstudio.
	
	There are two main directories:
	
		1) /data
		  - Houses raw data and/or summary data that underlies the study. 
			Subfolders correspond to different experimental goals (and
			therefore stimuli), and include /direction (direction tuning),
			/spatemp (spatiotemporal tuning). 
		  - For guidance on file naming 
			schemes, please see our scripts (next section) or send specific
			questions to VBB.
			
		2) /analyses_and_figures
		  - Code to recreate analyses and all figures located herein. It is
			strongly recommended to open `master_script.R` first to get 
			oriented and ensure that all necessary packages and generics are
			loaded. 
		  - Indiviudal scripts for sets of panels within each figure
			can be found within the corresponding figure's folder. Code for all
			anlyses can also be found in script for the corresponding figure
			panel(s). 
		  - Within a figure, panels found in a later position might
			rely on code for earlier panels (e.g., Fig 1E relies on analyses 
			for Fig 1B-D). It is therefore recommended to recreate panels in
			alphabetical order.
		  - Insertions of illustrations were largely performed via Adobe 
		    Illustrator. All illustrations can be found in the /illustrations
			sub-directory. 
	
	All work herein is licensed via an MIT license. See LICENSE file in root
	directory for further details.
		