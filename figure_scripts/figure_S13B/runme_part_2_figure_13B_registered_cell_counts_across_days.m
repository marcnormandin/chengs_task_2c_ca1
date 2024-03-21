% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This version makes Figure S13B, which are the cell counts
% and their identities across days (using registered cells). It tracks
% which deciles the cells change to.
%

close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER;
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    error('Run the script that makes the Excels first.');
end

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calData = load_data(calciumAnalysisResultsFolder);
calInfo = readtable(fullfile(LOCAL_OUTPUT_FOLDER, 'information_content_calcium_data.xlsx'));

%%
BA = calData.BestAligned;
CR = calData.CellRegTable;
registeredCellNames = cell(size(BA,1),1);

% Add the registered cell name to each row of BA
for iRow = 1:size(BA,1)
    cellName = BA.cellName{iRow};
    indMatch = find(ismember(CR.cellName, cellName));
    if length(indMatch) == 1
        registeredCellNames{iRow} = CR.registeredCellName{indMatch};
    end
end
BA.registeredCellName = registeredCellNames;
BA % to display it

%%
% Same as used for the decile scatter plot
binEdges = [min(BA.bestCorrelation), 0.31, 0.43, 0.51, 0.58, 0.64, 0.69, 0.74, 0.8, 0.85, 1];
BA(isnan(BA.bestCorrelation),:) = [];
binIds = discretize(BA.bestCorrelation, binEdges);
BA.binId = binIds;
BA(~any(BA.groupId == [1,2,3], 2),:) = [];
R = [];
rk = 1;
TAll = [];
for BIN_ID_A = 1:10
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    pairs = nchoosek(groupLabels, 2);
    for iPair = 1:size(pairs,1)
        groupLabelA = pairs(iPair,1);
        groupLabelB = pairs(iPair,2);
                
        BABinA = BA(BA.binId == BIN_ID_A & ismember(BA.groupLabel, groupLabelA),:);
        BABinB = BA(ismember(BA.groupLabel, groupLabelB),:);
        registeredCellNames = intersect(BABinA.registeredCellName, BABinB.registeredCellName);
        
        % Need to keep the indices into both lists
        nRegisteredCellNames = length(registeredCellNames);
        % Make a table
        T = [];
        k = 1;
        for iReg = 1:nRegisteredCellNames
            rcn = registeredCellNames{iReg};
        
            indA = find(ismember(BABinA.registeredCellName, rcn));
            indB = find(ismember(BABinB.registeredCellName, rcn));
        
            bestA = BABinA.bestCorrelation(indA);
            bestB = BABinB.bestCorrelation(indB);
            
            binIdA = BABinA.binId(indA);
            binIdB = BABinB.binId(indB); 
        
            T(k).iPair = iPair;
            T(k).registeredCellName = rcn;
            T(k).bestCorrelationA = bestA;
            T(k).bestCorrelationB = bestB;
            T(k).bestCorrelationChange = bestB - bestA;
            T(k).groupLabelA = groupLabelA;
            T(k).groupLabelB = groupLabelB;
            T(k).BIN_ID_A = BIN_ID_A;
            T(k).binIdA = binIdA;
            T(k).binIdB = binIdB;
            T(k).binIdChange = binIdB - binIdA;
            k = k + 1;
        end
        T = struct2table(T);
        TAll = [TAll; T];
    end
end

%% Save to do stats.
writetable(TAll, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_S13B_alternate_cellreg_count_results.xlsx'));


%% This computes percents as we vary the threshold for FS/FI
R1 = [];
rk = 1;
for BIN_ID_A = 1:10
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    pairs = nchoosek(groupLabels, 2);
    for iPair = 1:size(pairs,1)
        groupLabelA = pairs(iPair,1);
        groupLabelB = pairs(iPair,2);

        % Now calculate the percents
        T = TAll(ismember(TAll.groupLabelA, groupLabelA) & ismember(TAll.groupLabelB, groupLabelB) & TAll.BIN_ID_A <= BIN_ID_A,:);

        nCells = size(T,1);
        nRemained = sum(T.binIdA <= BIN_ID_A & T.binIdB <= BIN_ID_A);
        nChanged = nCells - nRemained;
        
        pRemained = nRemained / nCells * 100;
        pChanged = nChanged / nCells * 100;

        R1(rk).BIN_ID_A = BIN_ID_A;
        R1(rk).groupLabelA = groupLabelA;
        R1(rk).groupLabelB = groupLabelB;
        R1(rk).iPair = iPair;
        R1(rk).nCells = nCells;
        R1(rk).nRemained = nRemained;
        R1(rk).nChanged = nChanged;
        R1(rk).pRemained = pRemained;
        R1(rk).pChanged = pChanged;
        rk = rk + 1;
    end % iPair
end % iBin
R1 = struct2table(R1);
writetable(R1, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_S13B_alternate_cellreg_count_results_FS_remained_FS.xlsx'));

R2 = [];
rk = 1;
for BIN_ID_A = 1:10
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    pairs = nchoosek(groupLabels, 2);
    for iPair = 1:size(pairs,1)
        groupLabelA = pairs(iPair,1);
        groupLabelB = pairs(iPair,2);

        % Now calculate the percents
        T = TAll(ismember(TAll.groupLabelA, groupLabelA) & ismember(TAll.groupLabelB, groupLabelB) & TAll.BIN_ID_A <= BIN_ID_A,:);

        nCells = size(T,1);
        nRemained = sum(T.binIdA <= BIN_ID_A & T.binIdB > BIN_ID_A);
        nChanged = nCells - nRemained;
        
        pRemained = nRemained / nCells * 100;
        pChanged = nChanged / nCells * 100;

        R2(rk).BIN_ID_A = BIN_ID_A;
        R2(rk).groupLabelA = groupLabelA;
        R2(rk).groupLabelB = groupLabelB;
        R2(rk).iPair = iPair;
        R2(rk).nCells = nCells;
        R2(rk).nRemained = nRemained;
        R2(rk).nChanged = nChanged;
        R2(rk).pRemained = pRemained;
        R2(rk).pChanged = pChanged;
        rk = rk + 1;
    end % iPair
end % iBin
R2 = struct2table(R2);
writetable(R2, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_S13B_alternate_cellreg_count_results_FS_changed_to_FI.xlsx'));



%% This computes percents as we vary the threshold for FS
R1 = [];
rk = 1;
for BIN_ID_A = 1:10
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    pairs = nchoosek(groupLabels, 2);
    for iPair = 1:size(pairs,1)
        groupLabelA = pairs(iPair,1);
        groupLabelB = pairs(iPair,2);

        % Now calculate the percents
        T = TAll(ismember(TAll.groupLabelA, groupLabelA) & ismember(TAll.groupLabelB, groupLabelB) & TAll.BIN_ID_A >= BIN_ID_A,:);

        nCells = size(T,1);
        nRemained = sum(T.binIdA >= BIN_ID_A & T.binIdB >= BIN_ID_A);
        nChanged = nCells - nRemained;
        
        pRemained = nRemained / nCells * 100;
        pChanged = nChanged / nCells * 100;

        R1(rk).BIN_ID_A = BIN_ID_A;
        R1(rk).groupLabelA = groupLabelA;
        R1(rk).groupLabelB = groupLabelB;
        R1(rk).iPair = iPair;
        R1(rk).nCells = nCells;
        R1(rk).nRemained = nRemained;
        R1(rk).nChanged = nChanged;
        R1(rk).pRemained = pRemained;
        R1(rk).pChanged = pChanged;
        rk = rk + 1;
    end % iPair
end % iBin
R1 = struct2table(R1);
writetable(R1, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_S13B_alternate_cellreg_count_results_FI_remained_FI.xlsx'));


R2 = [];
rk = 1;
for BIN_ID_A = 1:10
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    pairs = nchoosek(groupLabels, 2);
    for iPair = 1:size(pairs,1)
        groupLabelA = pairs(iPair,1);
        groupLabelB = pairs(iPair,2);

        % Now calculate the percents
        T = TAll(ismember(TAll.groupLabelA, groupLabelA) & ismember(TAll.groupLabelB, groupLabelB) & TAll.BIN_ID_A >= BIN_ID_A,:);

        nCells = size(T,1);
        nRemained = sum(T.binIdA <= BIN_ID_A & T.binIdB < BIN_ID_A);
        nChanged = nCells - nRemained;
        
        pRemained = nRemained / nCells * 100;
        pChanged = nChanged / nCells * 100;

        R2(rk).BIN_ID_A = BIN_ID_A;
        R2(rk).groupLabelA = groupLabelA;
        R2(rk).groupLabelB = groupLabelB;
        R2(rk).iPair = iPair;
        R2(rk).nCells = nCells;
        R2(rk).nRemained = nRemained;
        R2(rk).nChanged = nChanged;
        R2(rk).pRemained = pRemained;
        R2(rk).pChanged = pChanged;
        rk = rk + 1;
    end % iPair
end % iBin
R2 = struct2table(R2);
writetable(R2, fullfile(LOCAL_OUTPUT_FOLDER, 'figure_S13B_alternate_cellreg_count_results_FI_changed_to_FS.xlsx'));






%% bestCorrelationChange
close all

hFig = figure('position', get(0, 'screensize'));
for iPair = 1:size(pairs,1)
    TAllPair = TAll(TAll.iPair == iPair,:);
    ax(iPair) = subplot(1,size(pairs,1),iPair);
    
    boxplot(abs(TAllPair.bestCorrelationChange), TAllPair.BIN_ID_A); %, 'boxstyle', 'filled')

    title(sprintf('%s -> %s', TAllPair.groupLabelA{1}, TAllPair.groupLabelB{2}))
    xlabel('Decile')
    ylabel(sprintf('Absolute Correlation Difference\n(%s - %s)', TAllPair.groupLabelB{2}, TAllPair.groupLabelA{1}))
    set(gca, 'fontsize', 18, 'fontweight', 'bold');
    grid on
end
linkaxes(ax, 'xy')
mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S13B_abs_correlation_differences_deciles', {'png', 'svg'});


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
    [MapsData, CellRegTable] = load_maps_data(folder, analysisSettings);

    % Verify the integrity of the data because we have to index between the tables.
    [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned);
    if n_errors > 0
        error('Data compromised.\n');
    end

    % Store it in a struct to make it easier to manage
    s = [];
    s.folder = folder;
    s.analysisSettings = analysisSettings;
    s.analysisResults = analysisResults;
    s.MapsData = MapsData;
    s.BestAligned = BestAligned;
    s.CellRegTable = CellRegTable;
end

function [MapsData, CellRegTable] = load_maps_data(folder, analysisSettings)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    MapsData = tmp.analysisInput.MapsData;
    CellRegTable = tmp.analysisInput.CellRegTable;
    % We need to get the mapping from sessions to "days"
    eliminateUnmatched = false;
    MapsData = helper_add_group_labels_to_table(MapsData, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
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
