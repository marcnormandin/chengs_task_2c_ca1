% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This reproduces Figure 7B (average rate difference bar plot). The full
% statistics were done in R.

close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;

LOCAL_OUTPUT_FOLDER = './local_output'; % for Figure S5B
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');
dataset = load_data(tetrodesAnalysisResultsFolder);
if dataset.analysisSettings.IS_CALCIUM_DATA
    dataset.datasetName = 'calcium';
else
    dataset.datasetName = 'tetrodes';
end

%% This step computes the mean rate differences for within and across per cell
% The stability types still need to be appended to this data
RMF = dataset.analysisResults.RateMatricesFiltered;
meanRateDiffWithin = zeros(size(RMF,1),1);
meanRateDiffAcross = zeros(size(RMF,1),1);
for iRow = 1:size(RMF,1)
    cellData = table2struct(RMF(iRow,:));
    rm = cellData.rateMatrix;
    [I,J] = meshgrid(cellData.rateMatrixContextIds, cellData.rateMatrixContextIds);
    
    [withinI, withinJ] = find(I == J);
    withinIJ = [withinI, withinJ];
    withinIJ(withinIJ(:,1) == withinIJ(:,2),:) = []; % eliminate the diagonal
    withinInds = sub2ind(size(I), withinIJ(:,1), withinIJ(:,2));
    
    [acrossI, acrossJ] = find(I ~= J);
    acrossIJ = [acrossI, acrossJ];
    acrossIJ(acrossIJ(:,1) == acrossIJ(:,2),:) = []; % eliminate the diagonal
    acrossInds = sub2ind(size(I), acrossIJ(:,1), acrossIJ(:,2));
    
    meanRateDiffWithin(iRow) = mean(rm(withinInds), 'all', 'omitnan');
    meanRateDiffAcross(iRow) = mean(rm(acrossInds), 'all', 'omitnan');
end
RMF.meanRateDiffWithin = meanRateDiffWithin;
RMF.meanRateDiffAcross = meanRateDiffAcross;

%% Append the labels
SessionsToGroups = load_sessions_to_groups_table(dataset.analysisSettings);
RMF = helper_add_group_labels_to_table(RMF, SessionsToGroups, true);

%% Format to save to Excel so that stats can be run in R
RMF_Save = RMF;
RMF_Save(:, ismember(RMF_Save.Properties.VariableNames, {'rateMatrix', 'rateMatrixTrialIds', 'rateMatrixContextIds'})) = []; % remove columns that are unneeded arrays
writetable(RMF_Save, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_7B_data.xlsx'));


%% Compute descriptive stats
[T] = compute_descriptive_stats(RMF);

% Save the descriptive stats (note, we did robust anova in R)
writetable(T, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_7B_descriptive_stats.xlsx'));

%% Make the figure
[hFig] = make_figure(T);
mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_7B_meanRateDifference', {'png', 'svg'});


%% Compute basic stats
function [T] = compute_descriptive_stats(RMF)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroupLabels = length(groupLabels);
    comparisonTypes = {'within', 'across'};
    nComparisonTypes = length(comparisonTypes);

    T = [];
    tk = 1;
    for iGroupLabel = 1:nGroupLabels
        groupLabel = groupLabels{iGroupLabel};
    
        for iComparisonType = 1:nComparisonTypes
            comparisonType = comparisonTypes{iComparisonType};
    
            x = get_data_subset(RMF, groupLabel, comparisonType);
            meanX = mean(x,'all','omitnan');
            stdX = std(x,1,'all','omitnan');
            nX = length(x);
            errX = stdX ./ sqrt(nX);

            T(tk).groupLabel = groupLabel;
            T(tk).comparisonType = comparisonType;
            T(tk).mean = meanX;
            T(tk).stderr = errX;
            T(tk).nCells = nX;
            tk = tk + 1;
        end
    end
    T = struct2table(T);
end % function


%% Make the bar plot figure
function [hFig] = make_figure(T)
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroupLabels = length(groupLabels);
    comparisonTypes = {'within', 'across'};
    nComparisonTypes = length(comparisonTypes);

    % Format the data for plotting
    meanX = nan(nComparisonTypes, nGroupLabels);
    errX = nan(size(meanX));
    for iGroupLabel = 1:nGroupLabels
        groupLabel = groupLabels{iGroupLabel};
        for iComparisonType = 1:nComparisonTypes
            comparisonType = comparisonTypes{iComparisonType};
            meanX(iComparisonType, iGroupLabel) = T.mean(ismember(T.groupLabel, groupLabel) & ismember(T.comparisonType, comparisonType));
            errX(iComparisonType, iGroupLabel) = T.stderr(ismember(T.groupLabel, groupLabel) & ismember(T.comparisonType, comparisonType));
        end
    end

    [colours] = colour_matrix_yellow_green();
    
    hFig = figure('position', get(0, 'screensize'));
    b = bar(1:3, meanX, 'FaceColor', 'flat');
    hold on
    for i = 1:size(meanX,1)
        b(i).CData = colours(i,:);
    end
    
    for i = 1:length(b)
       errorbar(b(i).XEndPoints, meanX(i,:), errX(i,:),'k.','linewidth',4);
    end
    grid on
    xticks(1:3);
    xticklabels(groupLabels);
    ylabel('Mean Firing Rate Difference [Hz]')
    legend(comparisonTypes, 'location', 'southoutside', 'orientation', 'horizontal');
    set(gca, 'fontsize', 18, 'fontweight', 'bold')
end % figure


%% Helpers
function [x] = get_data_subset(RMF, groupLabel, comparisonType)
    if strcmpi(comparisonType, 'within')
        varName = 'meanRateDiffWithin';
    elseif strcmpi(comparisonType, 'across')
        varName = 'meanRateDiffAcross';
    else
        error('unknown comparison type');
    end
    inds = find(ismember(RMF.groupLabel, groupLabel));
                   
    x = RMF.(varName)(inds);
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
    eliminateUnmatched = false;
    BestAligned = helper_add_group_labels_to_table(BestAligned, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function
