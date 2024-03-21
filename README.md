# This repository has MATLAB scripts that will reproduce most of the plots from our Cheng's Task 2C CA1 paper. The other repos reproduce the other figures.

The code requires the associated data files which are not stored in this repo.

How to set up:
- Download the MuLaNA repository
- Download CircStats2012 package
- Download the data repo associated with the paper
- Edit the contents of /figure_scripts/load_figure_config.m
- - Change the respective variables for the folders to point to the locations on your computer.

To reproduce the figures (from paper and supplement):
- Go to the associated figure name under the /figure_scripts folder
- Run the "runme_figure_*" script in that folder. Be sure to read the top of the file because some figures have multiple stages.
- The figures will be saved to the folder /figure_scripts/OUTPUT_FIGURES by default (can be changed in load_figure_config.m).

To re-compute some of the data
- Run /main_scripts/runme_calcium.m and /main_scripts/runme_tetrodes.m
- If you re-compute, then you will need to also update the figure configuration in /figure_scripts/load_figure_config.m
