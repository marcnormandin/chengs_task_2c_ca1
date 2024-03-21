% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Figure S12B and Figure S12C. This was written to perform a simple 
% per-cell per trial context prediction.

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

% To save the excel for further processing
LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

% These are the folders created by the main analysis script to store the
% results
calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');


calData = load_data(calciumAnalysisResultsFolder);
tetData = load_data(tetrodesAnalysisResultsFolder);

%% Settings that control the withheld algorithm
algoSettings = [];
algoSettings.limit_digs = {'Corr', 'Geo'};
algoSettings.min_maps = 2;

[calV] = predict_contexts_and_filter(calData, algoSettings);
calV.recording_type = cell(size(calV,1),1);
calV.recording_type = cellfun(@(x)('calcium'), calV.recording_type, 'UniformOutput',false);

[tetV] = predict_contexts_and_filter(tetData, algoSettings);
tetV.recording_type = cell(size(tetV,1),1);
tetV.recording_type = cellfun(@(x)('tetrode'), tetV.recording_type, 'UniformOutput',false);

% Combine
allV = [calV; tetV];

%% Export
cNames = {'animalName', 'sessionName', 'dayNum', 'cellName', 'animalSessionCellName', 'cellId', 'bestCorrelation', 'bestCorrelation_recomputed', 'isStable', 'groupId', 'groupLabel', 'prediction_accuracy', 'predict_n_total', 'predict_n_correct', 'isValid', 'recording_type'};
E = allV(:, ismember(allV.Properties.VariableNames, cNames));
outputFn = fullfile(LOCAL_OUTPUT_FOLDER, 'predict_contexts_using_withheld_map.xlsx');
writetable(E, outputFn)
fprintf('Results saved to %s\n', outputFn);

%%
function [V] = predict_contexts_and_filter(data, algoSettings)
    [T] = predict_contexts(data.MapsData, data.BestAligned, algoSettings.min_maps, algoSettings.limit_digs);
    
    % Filter out errors from T
    valid_rows = find(cellfun(@(x)isempty(x)==false, T.groupLabel));
    V = T(valid_rows,:);
    bad_rows = find(isnan(V.bestCorrelation_recomputed) | isnan(V.prediction_accuracy));
    V(bad_rows,:) = [];
    
    V(V.isValid == false,:) = [];
end

function [T] = predict_contexts(MapsData, BestAligned, min_maps, limit_digs)
    % Map a copy of BA so we don't have to resave the same data
    T = BestAligned; 
    T.bestCorrelation_recomputed = zeros(size(T,1),1);
    T.prediction_accuracy = nan(size(T,1),1);
    T.prediction_n_total = nan(size(T,1),1);
    T.prediction_n_correct = nan(size(T,1),1);
    T.prediction_true = cell(size(T,1),1);
    T.prediction_predicted = cell(size(T,1),1);
    T.isValid = false(size(T,1),1);
    T.averageMapContext1_recomputed = cell(size(T,1),1);
    T.averageMapContext2_recomputed = cell(size(T,1),1);
    T.mapsContext1_recomputed = cell(size(T,1),1);
    T.mapsContext2_recomputed = cell(size(T,1),1);


    for bestAligned_row_index = 1:size(BestAligned,1)
        animalName = BestAligned.animalName{bestAligned_row_index};
        sessionName = BestAligned.sessionName{bestAligned_row_index};
        cellName = BestAligned.cellName{bestAligned_row_index};
        bestAligned_cellDataInds = BestAligned.cellDataInds{bestAligned_row_index};
        cellData = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        mapsData_cellDataInds = cellData.cellDataInds;
        n_inds = length(bestAligned_cellDataInds);
        if n_inds ~= length(intersect(bestAligned_cellDataInds, mapsData_cellDataInds))
            continue
        end

        cellData.maps = cellData.maps ./ sum(cellData.maps, [1,2]);
        
        % Align the maps based on the best found rotation sequence
        amaps1 = align_maps(cellData.maps(:,:, cellData.contextIds == 1), BestAligned.rotationSequence1{bestAligned_row_index});
        amaps2 = align_maps(cellData.maps(:,:, cellData.contextIds == 2), BestAligned.rotationSequence2{bestAligned_row_index});
        amaps = amaps1;
        for i = 1:size(amaps2,3)
            amaps(:,:,size(amaps1,3)+i) = amaps2(:,:,i);
        end
        true_contexts = [zeros(size(amaps1,3),1); ones(size(amaps2,3),1)];

        if ~isempty(limit_digs)
            digs1 = cellData.digs(cellData.contextIds == 1);
            digs2 = cellData.digs(cellData.contextIds == 2);
            digs = digs1;
            for i = 1:length(digs2)
                digs(length(digs1)+i) = digs2(i);
            end

            valid_digs = find(ismember(digs, limit_digs));
            % filter
            amaps = amaps(:,:,valid_digs);
            true_contexts = true_contexts(valid_digs);
        end

        % check that the maps are valid
        is_valid = false(length(true_contexts),1);
        for i = 1:length(is_valid)
            map = squeeze(amaps(:,:,i));
            if isempty(map)
                continue;
            end
            if any(~isfinite(map), 'all')
                continue;
            end
            if all(map == 0, 'all')
                continue;
            end
            is_valid(i) = true;
        end

        amaps = amaps(:,:,is_valid);
        true_contexts = true_contexts(is_valid);
        if isempty(true_contexts)
            continue;
        end

        n1 = sum(true_contexts == 0);
        n2 = sum(true_contexts == 1);
        if n1 < min_maps || n2 < min_maps
            continue; % invalid
        end

        predicted_contexts = nan(size(true_contexts));
        n_trials = length(true_contexts);
        
        inds_all = 1:n_trials;
        for iTrial = 1:n_trials
            wmap = squeeze(amaps(:,:,iTrial));
            inds_not_withheld = setdiff(inds_all, iTrial);
            oamaps = amaps(:,:, inds_not_withheld);
            ocontexts = true_contexts(inds_not_withheld);
            c1wavg = mean(oamaps(:,:, ocontexts==0), 3);
            c2wavg = mean(oamaps(:,:, ocontexts==1), 3);
        
            c1corr = corrcoef(wmap(:), c1wavg);
            c1corr = c1corr(1,2);
        
            c2corr = corrcoef(wmap(:), c2wavg);
            c2corr = c2corr(1,2);
        
            if c1corr >= c2corr
                predicted_contexts(iTrial) = 0;
            else
                predicted_contexts(iTrial) = 1;
            end
        end % iTrial
        
        % Recomputed the best aligned correlation
        avgMap1 = mean(amaps(:,:, true_contexts==0),3);
        avgMap2 = mean(amaps(:,:, true_contexts==1),3);
        bestCorrelation_recomputed = corr(avgMap1(:), avgMap2(:));
        %BestAligned.bestCorrelation(bestAligned_row_index)
    
        % store
        T.bestCorrelation_recomputed(bestAligned_row_index) = bestCorrelation_recomputed;
        T.prediction_accuracy(bestAligned_row_index) = sum(true_contexts == predicted_contexts) / length(true_contexts);
        T.prediction_n_total(bestAligned_row_index) = length(true_contexts);
        T.prediction_n_correct(bestAligned_row_index) = sum(true_contexts == predicted_contexts);
        T.prediction_true{bestAligned_row_index} = true_contexts;
        T.prediction_predicted{bestAligned_row_index} = predicted_contexts;
        T.isValid(bestAligned_row_index) = true;
        T.averageMapContext1_recomputed{bestAligned_row_index} = avgMap1;
        T.averageMapContext2_recomputed{bestAligned_row_index} = avgMap2;

        % Store the maps as they were used for the average and analysis
        T.mapsContext1_recomputed{bestAligned_row_index} = amaps(:,:, true_contexts == 0);
        T.mapsContext2_recomputed{bestAligned_row_index} = amaps(:,:, true_contexts == 1);
    end % bestAligned_row_index
end % function



function [maps] = align_maps(maps, rotation_sequence)
    for i = 1:length(rotation_sequence)
        if rotation_sequence(i)==1
            m = maps(:,:,i);
            mrot = rot90(m,2);
            maps(:,:,i) = mrot;
            %fprintf('Rotated Map %d\n', i)
        end
    end
end


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
    [BestAligned, analysisSettings] = load_best_aligned(folder);
    [MapsData] = load_maps_data(folder, analysisSettings);

    % Verify the integrity of the data because we have to index between the tables.
    [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned);
    if n_errors > 0
        error('Data compromised.\n');
    end

    % Store it in a struct to make it easier to manage
    s = [];
    s.folder = folder;
    s.analysisSettings = analysisSettings;
    s.MapsData = MapsData;
    s.BestAligned = BestAligned;
end

function [MapsData] = load_maps_data(folder, analysisSettings)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    MapsData = tmp.analysisInput.MapsData;
    eliminateUnmatched = false;
    MapsData = helper_add_group_labels_to_table(MapsData, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function

function [BestAligned, analysisSettings] = load_best_aligned(folder)
    tmp = load(fullfile(folder, 'analysis_results.mat'));
    analysisSettings = tmp.analysisSettings;
    BestAligned = tmp.analysisResults.BestAligned;
    eliminateUnmatched = false;
    BestAligned = helper_add_group_labels_to_table(BestAligned, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function
