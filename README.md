## About this repository

This repository has MATLAB scripts that will produce plots for our Cheng's Task 2C CA1 neuroscience paper titled ***"Distinct neural mechanisms for heading retrieval and context recognition in the hippocampus during spatial reorientation"***. Other repos referenced in the paper will produce the figures not found here.

## Required Dataset files

The code requires the associated data files which are not stored in this repo. See the paper for a link to them.

## Required codes:
* Download MuLaNA (frozen for the paper): [https://github.com/marcnormandin/chengs_task_2c_ca1_MuLaNA](https://github.com/marcnormandin/chengs_task_2c_ca1_MuLaNA)
* Download Distinguishable Colors [https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors](https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors)
* Download Circular Statistics Toolbox: https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics
* Download boxplotGroupV2 [https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup](https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup)

## How to adapt the code for your system:

* Open the file `/figure_scripts/load_figure_config.m` in an editor.
* Change the respective folder variable values to be the locations on your local computer.

## To reproduce the figures (from paper and supplement):

* The code for each figure is located in subfolders of `/figure_scripts`. Navigate to the folder of the figure you wish to make and then run the associated main script; generally the main script will be prefixed `runme_figure`, but some scripts have have multiple stages so run those in sequence.

* The figures will be saved to the folder specified by `CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER`.
