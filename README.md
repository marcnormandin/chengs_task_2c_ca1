# This repository has MATLAB scripts that will produce plots for our Cheng's Task 2C CA1 paper. The other repos reproduce the other figures.

The code requires the associated data files which are not stored in this repo.

How to set up:
- Download MuLaNA: [https://github.com/marcnormandin/MuLaNA](https://github.com/marcnormandin/chengs_task_2c_ca1_MuLaNA/tree/main)
- Download Circular Statistics Toolbox: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
- Download the data repo associated with the paper
- Edit the contents of /figure_scripts/load_figure_config.m
- - Change the respective variables for the folders to point to the locations on your computer.

To reproduce the figures (from paper and supplement):
- Go to the associated figure name under the /figure_scripts folder
- Run the "runme_figure_*" script in that folder. Be sure to read the top of the file because some figures have multiple stages.
- The figures will be saved to the folder specified by CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER.
