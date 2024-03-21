% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Simple way to make the scripts run for individuals since the data will
% not be stored in the git repository and we need other repos. 

% All configuration variables are part of this structure.
CHENGS_TASK_2C_FIGURES_CONFIG = [];

% Set this to the folder where the data is located
CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER = 'D:\chengs_task_2c_ca1_project\data';

% Set this to the folder where the figures will be created.
CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER = 'D:\chengs_task_2c_ca1_project\output';

% Set these to where the github repos are saved
CHENGS_TASK_2C_FIGURES_CONFIG.MULANA_FOLDER = 'D:\chengs_task_2c_ca1_project\codes\MuLaNA';
CHENGS_TASK_2C_FIGURES_CONFIG.CIRC_STAT_FOLDER = 'D:\chengs_task_2c_ca1_project\codes\CircStat2012a';
CHENGS_TASK_2C_FIGURES_CONFIG.DISTINGUISHABLE_COLORS_FOLDER = 'D:\chengs_task_2c_ca_projects\codes\DistinguishableColors';
CHENGS_TASK_2C_FIGURES_CONFIG.MAIN_SCRIPTS_FOLDER = 'D:\chengs_task_2c_ca1_project\codes\chengs_task_2c_ca1';
CHENGS_TASK_2C_FIGURES_CONFIG.CIRCULAR_STATISTICS_TOOLBOX_FOLDER = 'D:\chengs_task_2c_ca1_project\codes\CircularStatisticsToolbox';
CHENGS_TASK_2C_FIGURES_CONFIG.BOXPLOTGROUP_FOLDER = 'D:\chengs_task_2c_ca1_project\codes\BoxPlotGroup';

% Add the paths we need. These have the common scripts we
% need in addition to the main analysis.
addpath(genpath(fullfile(CHENGS_TASK_2C_FIGURES_CONFIG.MAIN_SCRIPTS_FOLDER, 'main_scripts')));
addpath(genpath(CHENGS_TASK_2C_FIGURES_CONFIG.MULANA_FOLDER));
addpath(genpath(CHENGS_TASK_2C_FIGURES_CONFIG.CIRC_STAT_FOLDER));
addpath(genpath(CHENGS_TASK_2C_FIGURES_CONFIG.DISTINGUISHABLE_COLORS_FOLDER));
addpath(genpath(CHENGS_TASK_2C_FIGURES_CONFIG.CIRCULAR_STATISTICS_TOOLBOX_FOLDER));
addpath(genpath(CHENGS_TASK_2C_FIGURES_CONFIG.BOXPLOTGROUP_FOLDER));
