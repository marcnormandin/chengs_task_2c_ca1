% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This is the main script (for tetrodes) to compute the data tables we need to perform the
% further analysis.

close all
clear all
clc

run('../figure_scripts/load_figure_config.m')


addpath(genpath('.'));

tstart = tic;

dateTag = ml_util_gen_datetag();
OUTPUT_RESULTS_FOLDER = fullfile(pwd, 'analysis_results', sprintf('tetrodes_%s', dateTag));
if ~exist(OUTPUT_RESULTS_FOLDER, 'dir')
    mkdir(OUTPUT_RESULTS_FOLDER);
end




%% Configuration settings.
analysisSettings = analysis_settings_load();
analysisSettings.dateTag = dateTag;
analysisSettings.OUTPUT_RESULTS_FOLDER = OUTPUT_RESULTS_FOLDER;



%% Load the input. Mainly the rectangular and square maps.
analysisInput = analysis_load_celia_tetrode_input_data(analysisSettings);
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


function [analysisInput] = analysis_load_celia_tetrode_input_data(analysisSettings)
    % Load the main rectangular placemap data
    fprintf('Loading rectangular map data: ');
    [dataFilename, filePath] = uigetfile();
    tmp = load(fullfile(filePath, dataFilename));
    MapsData = load_celia_data_structure_as_table_v2(tmp.ratemaps.data);
    fprintf('done!\n');



    % Compute the square maps, which we need for the BFO 90
    fprintf('Computing square map data: ');
    sqMaps = mkSqRatemaps(tmp.ratemaps);
    MapsSquareData = load_celia_data_structure_as_table_v2(sqMaps);
    clear sqMaps
    clear tmp
    fprintf('done!\n');
        
    % If set, eliminate data where animal did not dig
    if analysisSettings.ELIMINATE_NO_DIGS
       inds = find(~ismember(MapsData.dig, {'Corr', 'Geo', 'Feat', 'Wrong'}));
       MapsData(inds,:) = [];

       inds = find(~ismember(MapsSquareData.dig, {'Corr', 'Geo', 'Feat', 'Wrong'}));
       MapsSquareData(inds,:) = [];
    end
    clear inds
    
    
    % INFORMATION RATE FILTER
    % Filter, if non-empty, using the information content
    if ~isempty(analysisSettings.TETRODES_INFORATE_FILTER_PERCENTILE)
        fprintf('Filtering out all maps whose information content (infoRate) is lower than the %d-th percentile.\n', analysisSettings.TETRODES_INFORATE_FILTER_PERCENTILE);
        % Rectangular maps
        infoRate = real(MapsData.infoRate(:));
        infoRateV = prctile(infoRate, analysisSettings.TETRODES_INFORATE_FILTER_PERCENTILE);
        
        MapsData(infoRate < infoRateV,:) = [];
        
        % We will filter the square maps by name because the square maps dont actually
        % have an infoRate
        uniqueAnimalSessionCellNames = unique(MapsData.animalSessionCellName);
        
        % Square maps
        badInds = find(~ismember(MapsSquareData.animalSessionCellName, uniqueAnimalSessionCellNames));
        MapsSquareData(badInds,:) = [];
    else
        fprintf('No filtering based on an infoRate threshold was applied.\n');
    end
    
    
    
    
    % Store results
    analysisInput = [];
    analysisInput.MapsData = MapsData;
    analysisInput.MapsSquareData = MapsSquareData;
    analysisInput.dataFilename = dataFilename;
end

