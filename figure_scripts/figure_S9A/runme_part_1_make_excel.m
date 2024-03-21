% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Excel files needed for part 2 (Figure S9) (but used Excel + R for paper).
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.DATA_FOLDER; %'../figure_S13B/local_output';
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

% To save the excel for further processing
LOCAL_OUTPUT_FOLDER = './local_output';
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

load(fullfile(INPUT_FOLDER,'ratemaps-2021_12_1_12_34_49_K1_d4.mat'))
TRedo = load_celia_data_structure_as_table_vIowa(ratemaps.data, true);
clear ratemaps

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');

calData = load_data(calciumAnalysisResultsFolder);
tetData = load_data(tetrodesAnalysisResultsFolder);

%%
MapsData = tetData.MapsData;
MapsData.infoRate_old = nan(size(MapsData,1),1);
MapsData.infoSpike_old = nan(size(MapsData,1),1);

for iRow = 1:size(MapsData,1)
    fprintf('Processing %d of %d\n', iRow, size(MapsData,1));
    animalSessionCellName = MapsData.animalSessionCellName{iRow};
    trialId = MapsData.trialId(iRow);
    % Find the match in the recalculated table
    ind = find(ismember(TRedo.animalSessionCellName, animalSessionCellName) & TRedo.trialId==trialId);
    if length(ind) ~= 1
        error('Something is wrong. Could not find %s in the recalculated table.', animalSessionCellName);
    end
    MapsData.mfr(iRow) = TRedo.mfr(ind);
    MapsData.pfr(iRow) = TRedo.pfr(ind);
    MapsData.infoRate(iRow) = TRedo.infoRate(ind);
    MapsData.infoSpike(iRow) = TRedo.infoSpike(ind);

    MapsData.infoRate_old(iRow) = TRedo.infoRate_old(ind);
    MapsData.infoSpike_old(iRow) = TRedo.infoSpike_old(ind);
end
tetData.MapsData = MapsData;
clear MapsData

%%
clc
close all
datasets = {calData, tetData};
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};
    [T] = compute(dataset);
    dataset.T = T;
    datasets{iDataset} = dataset;
end


%% Set the data as Isabel wants
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};
    T = dataset.T;
    E = T(:, ismember(T.Properties.VariableNames, {'animalName', 'sessionName', 'animalSessionCellName', 'isStable', 'groupId', 'groupLabel', 'bestCorrelation', 'ics_c1_avg', 'ics_c2_avg', 'ics_avg'}));
    E =  sortrows(E, {'isStable', 'groupLabel', 'animalName'});

    if dataset.analysisSettings.IS_CALCIUM_DATA
        outputFn = fullfile(LOCAL_OUTPUT_FOLDER, 'information_content_calcium_data.xlsx');
    else
        outputFn = fullfile(LOCAL_OUTPUT_FOLDER, 'information_content_tetrodes_data.xlsx');
    end
    writetable(E, outputFn);
    fprintf('Data written to %s\n', outputFn);
end





%% FUNCTIONS
function [T] = compute(dataset)
    [T] = compute_avg_ics_per_cell(dataset.analysisSettings, dataset.MapsData, dataset.BestAligned);
    SessionsToGroups = load_sessions_to_groups_table(dataset.analysisSettings);
    [T] = helper_add_group_labels_to_table(T, SessionsToGroups, true);
    T(~ismember(T.groupId, [1,2,3]),:) = [];
end % function




%% functions

function [T] = compute_avg_ics_per_cell(analysisSettings, MapsData, BestAligned)
    % Map a copy of BA so we don't have to resave the same data
    T = BestAligned; 
    % Just add the average ics to the best aligned table
    T.ics_c1_avg = nan(size(T,1),1);
    T.ics_c2_avg = nan(size(T,1),1);
    T.ics_avg = nan(size(T,1),1);

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

        c1inds = mapsData_cellDataInds(cellData.contextIds == 1);
        c2inds = mapsData_cellDataInds(cellData.contextIds == 2);

        % avg information per context
        if analysisSettings.IS_CALCIUM_DATA
            ic1_avg = mean(MapsData.ics(c1inds), 'omitnan');
            ic2_avg = mean(MapsData.ics(c2inds), 'omitnan');
            ic_avg = mean(MapsData.ics(mapsData_cellDataInds), 'omitnan');
        else % or infoSpike
            ic1_avg = mean(MapsData.infoRate(c1inds), 'omitnan');
            ic2_avg = mean(MapsData.infoRate(c2inds), 'omitnan');
            ic_avg = mean(MapsData.infoRate(mapsData_cellDataInds), 'omitnan');
        end


        % store
        T.ics_c1_avg(bestAligned_row_index) = ic1_avg;
        T.ics_c2_avg(bestAligned_row_index) = ic2_avg;
        T.ics_avg(bestAligned_row_index) = ic_avg;

        
    end % bestAligned_row_index
end % function

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
    [analysisResults, analysisSettings] = load_analysis_results(folder);
    tmp = load(fullfile(folder, 'analysis_results_group_averages.mat'));
    analysisResultsGroupAverages = tmp.analysisResultsGroupAverages;

    [BestAligned, analysisSettings2] = load_best_aligned(folder);
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
    s.analysisResults = analysisResults;
    s.analysisResultsGroupAverages = analysisResultsGroupAverages;
    s.MapsData = MapsData;
    s.BestAligned = BestAligned;
end

function [analysisResults, analysisSettings] = load_analysis_results(folder)
    tmp = load(fullfile(folder, 'analysis_results.mat'));
    analysisResults = tmp.analysisResults;
    analysisSettings = tmp.analysisSettings;
end % function



function [MapsData] = load_maps_data(folder, analysisSettings)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    MapsData = tmp.analysisInput.MapsData;
    fprintf('Loading the rectangular ratemaps.\n');
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


% This loads celias tetrode data into a form that can be used by the
% population vector code so the code can be general.
function [T] = load_celia_data_structure_as_table_vIowa(data, RECALC_INFORMATION)
    animalNames = extractfield(data, 'animal');
    numAnimals = length(animalNames);
    
    T = [];
    k = 1;
    for iAnimal = 1:numAnimals
        sessionsData = data(iAnimal).session;
        sessionNames = extractfield(sessionsData, 'name');
        animalName = data(iAnimal).animal;

        for iSession = 1:length(sessionNames)
            sessionName = sessionNames{iSession};
            trialData = sessionsData(iSession).trial;

            numTrials = length(trialData);

            % Get all of the unique cell names that are present for the
            % SESSION so we can give each a unique integer id
            cellNames = [];
            for iTrial = 1:numTrials
                tdata = trialData(iTrial);
                for iCell = 1:length(tdata.label)
                   %cellNames = cat(1, cellNames, tdata.label{iCell});
                   kk = length(cellNames)+1;
                   cellNames{kk} = tdata.label{iCell};
                end
            end
            uniqueCellNames = unique(cellNames);
            cellInds = arrayfun(@(x)(find(ismember(cellNames, uniqueCellNames{x}))), 1:length(uniqueCellNames), 'UniformOutput', false);
            uniqueCellIds = nan(length(cellNames), 1);
            for iCell = 1:length(cellInds)
                uniqueCellIds( cellInds{iCell} ) = iCell;
            end
            
            cnk = 1;
            for iTrial = 1:numTrials
                tdata = trialData(iTrial);

                dwellTimeMap = tdata.smoothedTime;
                probMap = dwellTimeMap ./ sum(dwellTimeMap, 'all', 'omitnan' );
                
                
                tid = tdata.trialNum;
                cid = tdata.context;
                numCells = length(tdata.label);
                for iCell = 1:numCells
                   T(k).animalName = animalName;
                   T(k).sessionName = sessionName;
                   T(k).dayNum = iSession;
                   T(k).cellName = tdata.label{iCell};
                   T(k).cellId = uniqueCellIds(cnk);
                   T(k).trialId = tid;
                   T(k).contextId = cid;

                   rateMap = tdata.maps{iCell};

                   T(k).map = rateMap;
                   
                   T(k).dig = tdata.dig;
                   
                   T(k).animalSessionCellName = sprintf('%s_%s_%s', T(k).animalName, T(k).sessionName, T(k).cellName);

                   if RECALC_INFORMATION
                       [informationRate, informationPerSpike] = ml_placefield_informationcontent_v2(rateMap, probMap);
                       T(k).mfr_old = tdata.mfr{iCell};
                       T(k).mfr = sum(probMap .* rateMap, 'all', 'omitnan');

                       T(k).pfr_old = tdata.pfr{iCell};
                       T(k).pfr = max(rateMap, [], 'all', 'omitnan');

                       T(k).infoRate_old = tdata.infoRate{iCell};
                       T(k).infoRate = informationRate;

                       T(k).infoSpike_old = tdata.infoSpike{iCell};
                       T(k).infoSpike = informationPerSpike;
                   else                  
                       % The square maps structure doesn't have all of the same
                       % field as the normal ratemaps structure for rectangular
                       % maps so we have to check.
                       if isfield(tdata, 'mfr')
                        T(k).mfr = tdata.mfr{iCell};
                       end
                       
                       if isfield(tdata, 'pfr')
                        T(k).pfr = tdata.pfr{iCell};
                       end
                       
                       if isfield(tdata, 'infoRate')
                        T(k).infoRate = tdata.infoRate{iCell};
                       end
                       
                       if isfield(tdata, 'infoSpike')
                        T(k).infoSpike = tdata.infoSpike{iCell};
                       end
                   end
                   
                   k = k + 1;
                   cnk = cnk + 1;
                end
            end
        end
    end
    T = struct2table(T);
end % function
