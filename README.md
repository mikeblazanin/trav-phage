# Paper information

## Fight not flight: parasites drive the bacterial evolution of resistance, not avoidance

Michael Blazanin, Jeremy Moore, Sydney Olsen, and Michael Travisano

**Abstract**
In the face of ubiquitous threats from parasites, hosts often evolve strategies to resist infection or to altogether avoid contact with parasites. At the microbial scale, bacteria frequently encounter viral parasites, bacteriophages. While bacteria are known to utilize a number of strategies to resist infection by phages, and can physically navigate their environment using complex motility behaviors, it is unknown whether bacteria evolve avoidance of phages. In order to answer this question, we combined experimental evolution and mathematical modeling. Experimental evolution of the bacterium Pseudomonas fluorescens in environments with differing spatial distributions of the phage Phi2 revealed that the host bacteria evolved resistance depending on parasite distribution and infectivity, but did not evolve dispersal to avoid parasite infection. Simulations using parameterized mathematical models of bacterial growth and swimming motility showed that this is a general finding: while increased dispersal is adaptive in the absence of parasites, in the presence of parasites that fitness benefit disappears and resistance becomes adaptive, regardless of the spatial distribution of parasites. Together, these experiments suggest that parasites should rarely, if ever, drive the evolution of bacterial avoidance via dispersal.

bioRxiv: https://doi.org/10.1101/2023.04.29.538831 

# Data and Analysis

All of the data and most of the analysis scripts are retained in the github repository https://github.com/mikeblazanin/trav-phage. 

Cleanup_Data.R is an R script that converts the raw data files stored in Raw_Data/ into consistent format data files stored in Clean_Data/. 

Analyze_Data.R is an R script that analyzes and visualized the clean data files from Clean_Data/.

Modeling data is generated by Matlab scripts located in the repo at https://github.com/jeremymoore558/ks_phage. Copies of the modeling data outputs are stored in the Modeling/ subdirectories. Visualization of the modeling data is done by Visualize_ModelOutput.R.

#Setup
To re-run the code in this repo, clone the repo then open the associated R Project file. The `renv` package should automatically bootstrap itself, downloading and installing the appropriate version of `renv` into the project library. Once this has completed, run
```
renv::restore()
```

This will install all package dependencies in the version used for our analyses. Then, you can run `Analyze_Data.R` or `Visualize_ModelOutput.R` to generate figures.
