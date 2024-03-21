% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This is the main script (for calcium) to compute the data tables we need to perform the
% further analysis.

close all
clear all
clc

run('../figure_scripts/load_figure_config.m')


addpath(genpath('.'));

tstart = tic;

dateTag = ml_util_gen_datetag();
OUTPUT_RESULTS_FOLDER = fullfile(pwd, 'analysis_results', sprintf('calcium_%s', dateTag));
if ~exist(OUTPUT_RESULTS_FOLDER, 'dir')
    mkdir(OUTPUT_RESULTS_FOLDER);
end

%% Configuration settings.
analysisSettings = analysis_settings_load();
analysisSettings.dateTag = dateTag;

analysisSettings.OUTPUT_RESULTS_FOLDER = OUTPUT_RESULTS_FOLDER;

analysisSettings.IS_CALCIUM_DATA = true;
analysisSettings.CALCIUM_MAPNAME_TO_USE = 'ilseMapsSmoothedAfter';


%% Load the input. Mainly the rectangular and square maps.
fprintf('Loading the calcium map data: \n');
analysisInput = analysis_load_marc_calcium_input_data(analysisSettings, []);
fprintf('\tFinished loading the calcium map data.\n'); 

save(fullfile(OUTPUT_RESULTS_FOLDER, 'analysis_input.mat'), 'analysisSettings', 'analysisInput')



%% Computations

% Main tables (no averaging)
analysisResults = main_analysis_compute_tables(analysisSettings, analysisInput);
save(fullfile(OUTPUT_RESULTS_FOLDER, 'analysis_results.mat'), 'analysisSettings', 'analysisResults');

%% Group averages (across days)
analysisResultsGroupAverages = main_analysis_group_averages_compute(analysisSettings, analysisResults);
save(fullfile(OUTPUT_RESULTS_FOLDER, 'analysis_results_group_averages.mat'), 'analysisSettings', 'analysisResultsGroupAverages');


%% Finish
telapsed = toc(tstart);
fprintf('Program completed in %0.2f minutes.\n\n', telapsed/60);

