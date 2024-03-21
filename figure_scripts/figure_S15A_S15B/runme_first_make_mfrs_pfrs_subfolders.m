% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Run this code first.
% This code computes the mfr and pfr for different occupancy thresholds and
% then saves the results in folders needed by the other scripts.

close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;
LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');
tetData = load_data(tetrodesAnalysisResultsFolder);

% Load the data for fraction of arena visited (per trial)
occupiedFractions = load(fullfile(INPUT_FOLDER, 'trial_map_fraction_occupied.mat'));
FractionOccupied = occupiedFractions.TTetrodes;

MapsData = tetData.analysisInput.MapsData;
BestAligned = tetData.BestAligned;

%%
MIN_FRACTION_OCCUPIED_LIST = 0.5:0.1:0.7;
RATE_MATRIX_RATETYPES = {'mfrs', 'pfrs'};

for iFraction = 1:length(MIN_FRACTION_OCCUPIED_LIST)

    MIN_FRACTION_OCCUPIED = MIN_FRACTION_OCCUPIED_LIST(iFraction);

    for iRateType = 1:length(RATE_MATRIX_RATETYPES)
        RATE_MATRIX_RATETYPE = RATE_MATRIX_RATETYPES{iRateType}; % 'pfrs' or 'mfrs'

        OUTPUT_FOLDER = fullfile(LOCAL_OUTPUT_FOLDER, RATE_MATRIX_RATETYPE, sprintf('min_fraction_occupied_%0.1f', MIN_FRACTION_OCCUPIED));
        if ~exist(OUTPUT_FOLDER, 'dir')
            mkdir(OUTPUT_FOLDER);
        end
        
        
        % Compute all of the rate matrices
        % There will be one row per unique cell
        algoSettings = tetData.analysisSettings; % This should be fine to also use for the calcium
        algoSettings.RATE_MATRIX_RATETYPE = RATE_MATRIX_RATETYPE; % changed from mfrs
        algoSettings.OUTPUT_RESULTS_FOLDER = OUTPUT_FOLDER; % So we don't corrupt the new results with the old.
        
        
        % Filter the MapsData
        keepRow = true(size(MapsData,1),1);
        for iRow = 1:size(MapsData,1)
            animalName = MapsData.animalName{iRow};
            sessionName = MapsData.sessionName{iRow};
            trialId = MapsData.trialId(iRow);
        
            % Look it up
            matchInd = find( ismember(FractionOccupied.animalName, animalName) & ...
                ismember(FractionOccupied.sessionName, sessionName) & ...
                FractionOccupied.trialId == trialId );
            if length(matchInd) == 1
                fo = FractionOccupied.fractionOccupied(matchInd);
                if fo <= MIN_FRACTION_OCCUPIED
                    keepRow(iRow) = false;
                end
            else
                fprintf('Strange: %s %s %d has %d matches\n', animalName, sessionName, trialId, length(matchInd));
            end
        end
        MapsDataFiltered = MapsData(keepRow,:);
        
        [analysisResults, analysisResultsGroupAverages] = run_ratematrix_analysis(MapsDataFiltered, BestAligned, algoSettings);
        make_ratematrix_plots(analysisResultsGroupAverages, algoSettings)

        SessionsToGroups = load_sessions_to_groups_table(algoSettings);
        SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];
        analysisResults.RateMatricesFiltered = helper_add_group_labels_to_table(analysisResults.RateMatricesFiltered, SessionsToGroups, true);

        RM = analysisResults.RateMatricesFiltered;
        
        
        % Add the best aligned correlations associated with the matrices
        RM.bestCorrelation = nan(size(RM,1),1);
        for iRow = 1:size(RM,1)
            animalSessionCellName = RM.animalSessionCellName{iRow};
            match = BestAligned.bestCorrelation(ismember(BestAligned.animalSessionCellName, animalSessionCellName));
            if length(match) ~= 1
                error('Someting wrong');
            end
            RM.bestCorrelation(iRow) = match;
        end
        
        
        T = [];
        k = 1;
        for iRow = 1:size(RM,1)
            rm = RM.rateMatrix{iRow};
            rateMatrixContextIds = RM.rateMatrixContextIds(iRow,:);
            nTrials = length(rateMatrixContextIds);
        
            within_values = [];
            across_values = [];
            for iTrialA = 1:nTrials
                for iTrialB = 1:nTrials
                    contextA = rateMatrixContextIds(iTrialA);
                    contextB = rateMatrixContextIds(iTrialB);
                    v = rm(iTrialA, iTrialB);
                    if contextA == contextB
                        within_values = [within_values, v];
                    else
                        across_values = [across_values, v];
                    end
                end
            end
        
            T(k).animalName = RM.animalName{iRow};
            T(k).sessionName = RM.sessionName{iRow};
            T(k).dayLabel = RM.groupLabel{iRow};
            T(k).dayId = RM.groupId(iRow);
            T(k).cellName = RM.cellName{iRow};
            T(k).animalSessionCellName = RM.animalSessionCellName{iRow};
        
            T(k).RATE_NAME = RM.rateName{iRow};
            T(k).MIN_FRACTION_OCCUPIED = MIN_FRACTION_OCCUPIED;
            T(k).DATASET_NAME = 'tetrodes';
        
            T(k).bestCorrelation = RM.bestCorrelation(iRow);
            if RM.bestCorrelation(iRow) <= 0.3
                T(k).stability = 'unstable';
            else
                T(k).stability = 'stable';
            end
        
            T(k).within_mean = mean(within_values, 'all', 'omitnan');
            T(k).within_std = std(within_values, 1, 'all', 'omitnan');
            T(k).within_median = median(within_values, 'all', 'omitnan');
        
            T(k).across_mean = mean(across_values, 'all', 'omitnan');
            T(k).across_std = std(across_values, 1, 'all', 'omitnan');
            T(k).across_median = median(across_values, 'all', 'omitnan');
            
            k = k + 1;
        end % iRow
        T = struct2table(T);

        % Make all the plots (PER CELL)
        stabilityTypes = {'any', 'unstable', 'stable'};
        for iStabilityType = 1:length(stabilityTypes)
            stabilityType = stabilityTypes{iStabilityType};
    
            [xMean, xStd, xN, xErr, groupLabels, comparisonTypes] = average_T_percell(T, stabilityType);
            [hFig] = plot_T_precalc(xMean, xErr, stabilityType, groupLabels, comparisonTypes);
            mulana_savefig(hFig, OUTPUT_FOLDER, sprintf('fig_barplot_percell_%s', stabilityType), {'png', 'svg', 'fig'});
            close(hFig);
        end % iStabilityType


        S = sortrows(T, {'stability', 'animalName', 'dayLabel', 'sessionName', 'animalSessionCellName'});

        writetable(S, fullfile(OUTPUT_FOLDER, sprintf('tetrodes_rate_matrix_within_across_per_cell_%s_min_fraction_occupied_%0.1f.xlsx', RATE_MATRIX_RATETYPE, MIN_FRACTION_OCCUPIED)));
    end % MIN_FRACTION_OCCUPIED
end % RATE_MATRIX_RATETYPE

fprintf('End of program\n');

%%

function [hFig] = plot_T_precalc(xMean, xErr, stabilityType, groupLabels, comparisonTypes)
    hFig = figure('position', get(0, 'screensize'));
    b = bar(1:3, xMean);
    hold on
    for i = 1:length(b)
       errorbar(b(i).XEndPoints, xMean(i,:), xErr(i,:),'k.','linewidth',4);
    end
    %errorbar(1:3, meanX(1,:), errX(1,:));
    grid on
    xticks(1:3);
    xticklabels(groupLabels);
    ylabel('Mean Firing Rate Difference [Hz]')
    legend(comparisonTypes, 'location', 'southoutside', 'orientation', 'horizontal');
    title(sprintf('Tetrodes (%s) averaged per cell', stabilityType), 'interpreter', 'none');
    set(gca, 'fontsize', 18, 'fontweight', 'bold')
end % function


function [x] = get_data_from_T(T, groupLabel, comparisonType, stabilityType)
    if strcmpi(stabilityType,'any')
        inds = find(ismember(T.dayLabel, groupLabel));
    else
        inds = find(ismember(T.dayLabel, groupLabel) & ismember(T.stability, stabilityType));
    end
    
    varName = sprintf('%s_mean', comparisonType);
    
    x = T.(varName)(inds);
end

function [xMean, xStd, xN, xErr, groupLabels, comparisonTypes] = average_T_percell(T, stabilityType)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);
    comparisonTypes = {'within', 'across'};
    nComparisonTypes = length(comparisonTypes);
    xMean = zeros(nComparisonTypes, nGroups);
    xStd = zeros(size(xMean));
    xN = zeros(size(xMean));
    xErr = zeros(size(xMean));
    
    for iComparisonType = 1:nComparisonTypes
        comparisonType = comparisonTypes{iComparisonType};
    
        for iGroup = 1:nGroups
            groupLabel = groupLabels{iGroup};
    
            [x] = get_data_from_T(T, groupLabel, comparisonType, stabilityType);
    
            xMean(iComparisonType, iGroup) = mean(x, 'all', 'omitnan');
            xStd(iComparisonType, iGroup) = std(x,1,'all','omitnan');
            xN(iComparisonType, iGroup) = sum(isfinite(x),'all');
            xErr(iComparisonType, iGroup) = xStd(iComparisonType, iGroup) ./ sqrt(xN(iComparisonType, iGroup));
        end
    end
end % function




function make_ratematrix_plots(analysisResultsGroupAverages, algoSettings)
    % Plot
    plottingSettings = plotting_settings_load(algoSettings);
    plottingSettings.OUTPUT_FOLDER = algoSettings.OUTPUT_RESULTS_FOLDER;
    
    
    plot_and_save_rate_matrices_and_hists_averaged_days123(algoSettings, analysisResultsGroupAverages, plottingSettings);
    plot_and_save_stable_rate_matrices_and_hists_averaged_days123(algoSettings, analysisResultsGroupAverages, plottingSettings);
end

function [analysisResults, analysisResultsGroupAverages] = run_ratematrix_analysis(MapsData, BestAligned, algoSettings)

    analysisResults = [];

    SessionsToGroups = load_sessions_to_groups_table(algoSettings);
    SessionsToGroups(~ismember(SessionsToGroups.dayNum, [1,2,3]),:) = [];
    

    [analysisResults.RateMatrices, success] = compute_rate_matrix_table(MapsData, algoSettings.RATE_MATRIX_RATETYPE, algoSettings.RATE_MATRIX_NUM_TRIALS);
    if ~success
        fprintf('Error while processing Rate Matrices (%s)\n', algoSettings.RATE_MATRIX_RATETYPE);
    end
    
    % Apply the filtering to tables if applicable
    analysisResults.RateMatricesFiltered = analysisResults.RateMatrices(analysisResults.RateMatrices.numRates >= algoSettings.RATE_MATRIX_FILTER_MIN_TRIALS,:);
    
    % COMPUTE AVERAGES
    analysisResultsGroupAverages = [];
    RateMatricesAveraged = compute_per_group_animal_averaged_rate_difference_matrix(SessionsToGroups, analysisResults.RateMatricesFiltered, algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL);
    
    % STORE
    analysisResultsGroupAverages.RateMatrices = RateMatricesAveraged;
    
    % COMPUTE AND STORE
    analysisResultsGroupAverages.RateMatricesHistogram = compute_group_average_rate_differences_within_across_for_table(RateMatricesAveraged);
    
    
    % STABILITY rate matrices and histograms for days 1,2,3

    % Compute
    [RMAverageStable, RMAverageUnstable, RMAverageStableWithinAcross, RMAverageUnstableWithinAcross] = compute_stability_animal_averaged_rate_matrices_within_across(SessionsToGroups, analysisResults.RateMatricesFiltered, BestAligned, ...
        algoSettings.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL);
    
    % STORE
    analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStable = RMAverageStable;
    analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstable = RMAverageUnstable;
    analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStableWithinAcross = RMAverageStableWithinAcross;
    analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstableWithinAcross = RMAverageUnstableWithinAcross;
end % function

function plot_and_save_rate_matrices_and_hists_averaged_days123(algoSettings, analysisResultsGroupAverages, plottingSettings)
    if ~algoSettings.IS_CALCIUM_DATA
        % Plotting
        RateMatricesAveraged = analysisResultsGroupAverages.RateMatrices;
        % Plot only days 1,2,3
        RateMatricesAveraged(~ismember(RateMatricesAveraged.groupId, [1,2,3]),:) = [];
        
        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RateMatricesAveraged, plottingSettings.RATE_MATRIX_CLIM);
        
        fnPrefix = sprintf("%srate_matrices_animalaveraged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end

    % Rate bars
    if ~algoSettings.IS_CALCIUM_DATA
        RateMatricesAveragedHistogram = analysisResultsGroupAverages.RateMatricesHistogram;
        hFig = plot_average_rate_differences_within_across_bars_table(RateMatricesAveragedHistogram);
        title('ANIMAL AVERAGED');
        fnPrefix = sprintf("%srate_hist_animalaveraged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end

end % function


function plot_and_save_stable_rate_matrices_and_hists_averaged_days123(algoSettings, analysisResultsGroupAverages, plottingSettings)
    if ~algoSettings.IS_CALCIUM_DATA
        RMAverageStable = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStable;
        RMAverageUnstable = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstable;
        RMAverageStableWithinAcross = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageStableWithinAcross;
        RMAverageUnstableWithinAcross = analysisResultsGroupAverages.StabilityRateMatrices.RateMatrixAverageUnstableWithinAcross;

        % Plot
        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RMAverageStable, plottingSettings.RATE_MATRIX_CLIM);
        sgtitle(sprintf('Stable (Per animal (%s, %s))', algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));    
        fnPrefix = sprintf("%sstable_rate_matrices_animalaveraged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);


        hFig = plot_per_group_animal_averaged_rate_difference_matrix_table(RMAverageUnstable, plottingSettings.RATE_MATRIX_CLIM);
        sgtitle(sprintf('Unstable (Per animal (%s, %s))', algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sunstable_rate_matrices_animalaveraged_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);

        % Hist 1
        hFig = plot_average_rate_differences_within_across_bars_table(RMAverageStableWithinAcross);
        sgtitle(sprintf('Stable (Per animal (%s, %s))', algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sstable_rate_matrices_animalaveraged_hist_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);


        % Hist 2
        hFig = plot_average_rate_differences_within_across_bars_table(RMAverageUnstableWithinAcross);
        sgtitle(sprintf('Unstable (Per animal (%s, %s))', algoSettings.RATE_MATRIX_NORMALIZATION_PER_CELL, algoSettings.RATE_MATRIX_NORMALIZATION_PER_ANIMAL));
        fnPrefix = sprintf("%sunstable_rate_matrices_animalaveraged_hist_days123", plottingSettings.FIGURE_PREFIX);
        ml_savefig(hFig, plottingSettings.OUTPUT_FOLDER, fnPrefix, plottingSettings.FIGURE_FORMATS_TO_SAVE);
        close(hFig);
    end % not calcium
end % function

%% Loading functions
function [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned)
    n_errors = 0;
    bad_ba_rows = [];
    for bestAligned_row_index = 1:size(BestAligned,1)
        animalName = BestAligned.animalName{bestAligned_row_index};
        sessionName = BestAligned.sessionName{bestAligned_row_index};
        cellName = BestAligned.cellName{bestAligned_row_index};
        bestAligned_cellDataInds = BestAligned.cellDataInds{bestAligned_row_index};
        [output] = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        mapsData_cellDataInds = output.cellDataInds;
        n_inds = length(bestAligned_cellDataInds);
        if n_inds ~= length(intersect(bestAligned_cellDataInds, mapsData_cellDataInds))
            %fprintf('Error with best aligned row index: %d\n', bestAligned_row_index);
            n_errors = n_errors + 1;
            bad_ba_rows(n_errors) = bestAligned_row_index;
        end
    end % bestAligned_row_index

end


function [s] = load_data(folder)
    [analysisResults, BestAligned, analysisSettings] = load_best_aligned(folder);
    %[MapsData] = load_maps_data(folder, algoSettings);
    [analysisInput] = load_analysis_input(folder);

    % Verify the integrity of the data because we have to index between the tables.
    [n_errors, bad_ba_rows] = verify_data_integrity(analysisInput.MapsData, BestAligned);
    if n_errors > 0
        error('Data compromised.\n');
    end

    % Store it in a struct to make it easier to manage
    s = [];
    s.folder = folder;
    s.analysisSettings = analysisSettings;
    s.analysisResults = analysisResults;
    s.analysisInput = analysisInput;
    s.BestAligned = BestAligned;
end


function [analysisInput] = load_analysis_input(folder)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    analysisInput = tmp.analysisInput;
end % function


function [analysisResults, BestAligned, analysisSettings] = load_best_aligned(folder)
    tmp = load(fullfile(folder, 'analysis_results.mat'));
    analysisSettings = tmp.analysisSettings;
    analysisResults = tmp.analysisResults;
    BestAligned = tmp.analysisResults.BestAligned;
    % We need to get the mapping from sessions to "days"
    eliminateUnmatched = false;
    BestAligned = helper_add_group_labels_to_table(BestAligned, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function




