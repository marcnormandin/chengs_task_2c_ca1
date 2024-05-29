% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure S4, which is of the BFO90 results for the dig pairs.
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

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
tetrodesAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');

calData = load_data(calciumAnalysisResultsFolder);
tetData = load_data(tetrodesAnalysisResultsFolder);

datasets = {calData, tetData};
%clear calData tetData

%% Recompute the bf90 per cell so that we can keep track of the associated dig pairs that we need
clc
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};
    NormalPerCell90 = recompute_bfo90_percell(dataset);
    isCalciumData = dataset.analysisSettings.IS_CALCIUM_DATA;
    if isCalciumData
        NormalPerCell90.isCalciumData = true(size(NormalPerCell90,1),1);
    else
        NormalPerCell90.isCalciumData = false(size(NormalPerCell90,1),1);
    end
    dataset.NormalPerCell90 = NormalPerCell90;
    datasets{iDataset} = dataset;
end % iDataset

%% Create a single results table
clc
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};

    R = compute_angle_behaviour_matrices(dataset);

    dataset.R = R;
    datasets{iDataset} = dataset;
end % iDataset

%% Compute the averages for each day (for each dataset)
clc
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroups = length(groupLabels);
digTypes = {'Corr', 'Feat', 'Geo', 'Wrong'};
digTypesPlot = {'Corr', 'Near', 'Geo', 'Far'};
nDigTypes = length(digTypes);
rotations = [0, 90, 180, 270];
nRotations = length(rotations);


for iDataset = 1:length(datasets)
    F = [];
    k = 1;

    dataset = datasets{iDataset};

    R = dataset.R;

    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        RG = R(ismember(R.groupLabel, groupLabel),:);

        BG = zeros(nDigTypes, nDigTypes, nRotations);
        for iAnimal = 1:size(RG,1)
            BG = BG + RG.B_normed_dim12{iAnimal};
        end
        F(k).B = BG / size(RG,1); % 2023-08-9. Added dividing by number of animals
        F(k).rotations = RG.rotations(1,:);
        F(k).digTypes = digTypes;
        F(k).groupId = RG.groupId(1);
        F(k).groupLabel = groupLabel;
       
        k = k + 1;

    end % iGroup
    F = struct2table(F);
    dataset.F = F;

    datasets{iDataset} = dataset;
end % iDataset

%% Separate each 3d matrix (dataset.F) into a single row
clc
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};

    F = dataset.F;
    R = [];
    k = 1;

    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};

        imatch = find(ismember(F.groupLabel, groupLabel));
        B = F.B{imatch};

        % this is hacky, but I don't care
        rotations = F.rotations(imatch,:);
        S0 = squeeze(B(:,:,rotations==0));
        S90 = squeeze(B(:,:,rotations==90));
        S180 = squeeze(B(:,:,rotations==180));
        S270 = squeeze(B(:,:,rotations==270));

        % 0 and 180 is just a straight sum since the order doesn't change
        % the angle group
        T0 = make_matrix_triangular(S0);
        T180 = make_matrix_triangular(S180);

        % The 90 and 270 swap the lower triangle of values
        T90 = zeros(size(S90));
        T270 = zeros(size(S270));        
        for i = 1:size(T90,1) % size of the matrices are the same so just do both at the same time.
            for j = i:size(T90,2)
                if j > i % Lower half (below diagonal)
                    % 
                    T90(j,i) = S90(j,i) + S270(i,j);
                    T270(j,i) = S270(j,i) + S90(i,j);

                    T90(i,j) = nan;
                    T270(i,j) = nan;
                else
                    % Upper half + diagonal (equal digs)
                    T90(i,j) = S90(i,j);
                    T270(i,j) = S270(i,j);
                end
            end
        end

        % Store the results
        for iRotation = 1:nRotations
            M = squeeze(B(:,:, iRotation));

            if F.rotations(imatch,iRotation) == 0
                MT = T0;
            elseif F.rotations(imatch, iRotation) == 90
                MT = T90;
            elseif F.rotations(imatch, iRotation) == 180
                MT = T180;
            elseif F.rotations(imatch, iRotation) == 270
                MT = T270;
            else
                error('asdfadsfafadsfasfasfas')
            end
            R(k).rotation = F.rotations(imatch,iRotation);
            R(k).digTypes = F.digTypes(imatch,:);
            R(k).groupId = F.groupId(imatch);
            R(k).groupLabel = F.groupLabel{imatch};
            R(k).matrix = M;
            R(k).matrix_triangular = MT;
            R(k).matrix_normed_deg = M ./ sum(M, 'all', 'omitnan');
            k = k + 1;
        end
    end
    dataset.matrices = struct2table(R);

    % This in inefficient but I dont care
    matrices = dataset.matrices;
    matrices.matrix_triangular_prob = cell(size(matrices,1),1);
    for iRow = 1:size(matrices,1)

        groupLabel = matrices.groupLabel{iRow};
        
        % Get the total matrix (we should have 4 rows of results)
        imatches = find(ismember(matrices.groupLabel, groupLabel));
        totalMatrix = zeros(nDigTypes, nDigTypes);
        for iInd = 1:length(imatches)
            totalMatrix = totalMatrix + matrices.matrix_triangular{imatches(iInd)};
        end

        % 
        MT = matrices.matrix_triangular{iRow} ./ totalMatrix;
        matrices.matrix_triangular_prob{iRow} = MT;

        % Norm each triangular matrix within a deg
        MTT = matrices.matrix_triangular{iRow};
        matrices.matrix_triangular_normed_deg{iRow} =  MTT ./ sum(MTT, 'all', 'omitnan');
    end % iRow



    dataset.matrices = matrices;
    datasets{iDataset} = dataset;
end






%% Plot combined (tetrodes + calcium)  Figure S3
clc
T = [];
tk = 1;
% First just go through and make the combined matrices, then plot after
for iDataset = 1:length(datasets)
    dataset = datasets{iDataset};
    if dataset.analysisSettings.IS_CALCIUM_DATA
        datasetName = 'calcium';
    else
        datasetName = 'tetrodes';
    end

    F = dataset.F;
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};
        B = F.B{ismember(F.groupLabel, groupLabel)};
        for iRotation = 1:nRotations
            M = squeeze(B(:,:, iRotation));

            T(tk).datasetName = datasetName;
            T(tk).groupLabel = groupLabel;
            T(tk).rotationDeg = rotations(iRotation);
            T(tk).matrix = M;
            tk = tk + 1;
        end
    end % iGroup
    
end % iDataset
T = struct2table(T)

% Now average
TA = [];
tk = 1;
for iGroup = 1:nGroups
    groupLabel = groupLabels{iGroup};
    for iRotation = 1:nRotations
        rotationDeg = rotations(iRotation);
        % subset (2 matrices)
        S = T(ismember(T.groupLabel, groupLabel) & T.rotationDeg == rotationDeg,:);
        SM = zeros(nDigTypes, nDigTypes, size(S,1));
        for k = 1:size(S,1)
            SM(:,:,k) = S.matrix{k};
        end
        TA(tk).groupLabel = groupLabel;
        TA(tk).rotationDeg = rotationDeg;
        TA(tk).datasetAveragedMatrix = squeeze(mean(SM,3, 'omitnan'));
        tk = tk + 1;
    end
end
TA = struct2table(TA)


%% Make the Figure S3
hFig = figure('position', get(0, 'screensize'));
ax = [];
k = 1;
PP = nGroups;
QQ = nRotations;
clims = zeros(PP*QQ,2);
% digTypes2 is new names for the dig areas

for iGroup = 1:nGroups
    groupLabel = groupLabels{iGroup};

    for iRotation = 1:nRotations
        ax(k) = subplot(PP,QQ, (iGroup-1)*nRotations + iRotation);
        rotationDeg = rotations(iRotation);

        k = k + 1;
        M = TA.datasetAveragedMatrix{ismember(TA.groupLabel, groupLabel) & TA.rotationDeg == rotationDeg};

        P = nan(size(M,1)+1, size(M,2)+1);
        P(1:size(M,1), 1:size(M,2)) = M;
        pcolor(P);
        shading flat
        colormap parula
        xticks(1.5:1:4.5)
        xticklabels(digTypesPlot);
        yticks(1.5:1:4.5);
        yticklabels(digTypesPlot);
        xlabel('Second Dig')
        ylabel('First Dig')
        title(sprintf('%s: %d deg', groupLabel, rotations(iRotation)));
        colorbar

        cb=colorbar;
        t=get(cb,'Limits');
        set(cb,'Ticks',linspace(t(1),t(2),3))


        clims(k,:) = clim;
    end
end % iGroup
clim_min = 0;
clim_max = max(clims, [], 'all');

USE_SAME_CLIMS = false;
if USE_SAME_CLIMS
    % make the colours represent the same values
    for iGroup = 1:nGroups
        groupLabel = groupLabels{iGroup};

        for iRotation = 1:nRotations
            rotationDeg = rotations(iRotation);
            ax(k) = subplot(PP,QQ, (iGroup-1)*nRotations + iRotation);
            clim([clim_min, clim_max]);
        end
    end % iGroup
end

sgtitle('Averaged Calcium + Tetrode Plots');
mulana_savefig(hFig, OUTPUT_FOLDER, sprintf('figure_S4_averaged_bfo90_digtype_matrices'), {'svg', 'png', 'fig'})


%% Added to save excel for NatComms requirement
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
rotationDegs = [0, 90, 180, 270];
X = [];
k = 1;
for iRow = 1:size(TA,1)
    
    XM = TA.datasetAveragedMatrix{iRow};
    for i = 1:length(digTypesPlot)
        firstDig = digTypesPlot{i};
        for j = 1:length(digTypesPlot)
            secondDig = digTypesPlot{j};

            X(k).dayName = TA.groupLabel{iRow};
            X(k).rotationDeg = TA.rotationDeg(iRow);
            X(k).currentDig = firstDig;
            X(k).subsequentDig = secondDig;
            X(k).probability = XM(i,j);
            k = k + 1;
        end
    end
end
X = struct2table(X)
writetable(X, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_S4.xlsx'), 'Sheet', 'figure_S4');

%% FUNCTIONS
function [MT] = make_matrix_triangular(M)
    % This doesn't take into account that a 90 for one pair will be 270 the
    % other
    MT = M + M' - M.*eye(4);
    for i = 1:size(MT,1)
        for j = i+1:size(MT,2)
            MT(i,j) = nan;
        end
    end
end % function

% This is the main function that does the counts
function [R] = compute_angle_behaviour_matrices(dataset)
    R = []; % results table
    k = 1;
    NPC90 = dataset.NormalPerCell90;
    
    digTypes = unique(dataset.analysisInput.MapsSquareData.dig);
    if ~all(ismember(digTypes, {'Corr', 'Feat', 'Geo', 'Wrong'}))
        error('Dig types are not the ones we want.');
    end
    
    Sessions = get_sessions_from_table(NPC90);
    for iSession = 1:size(Sessions,1)
        animalName = Sessions.animalName{iSession};
        sessionName = Sessions.sessionName{iSession};
        % Get percell90 associated with cells of the selected session
        T = NPC90(ismember(NPC90.animalName, animalName) & ismember(NPC90.sessionName, sessionName),:);
    
        nDigTypes = length(digTypes);
        rotations = [0, 90, 180, 270];
        nRotations = length(rotations);
        
        B = nan(nDigTypes, nDigTypes, nRotations, size(T,1)); % 3d, each plane is for one angle
        
        contextType = 'any';
        BDigs = nan(nDigTypes, nDigTypes, size(T,1)); % 3rd is cell index
        for iRow = 1:size(T,1) % for every cell of the session
            digpairs = T.(sprintf('digpairs_%s', contextType)){iRow};
            rotationInds = T.(sprintf('vind_%s', contextType)){iRow};
    
            for kk = 1:size(digpairs,1)
                digA = digpairs{kk,1};
                digAInd = find(ismember(digTypes, digA));
    
                digB = digpairs{kk,2};
                digBInd = find(ismember(digTypes, digB));
    
                rotInd = rotationInds(kk);
    
                % add to counts
                prev = B(digAInd, digBInd, rotInd, iRow);
                if isnan(prev)
                    prev = 0;
                end
                B(digAInd, digBInd, rotInd, iRow) = prev + 1;

                % add for the total of digs
                prev = BDigs(digAInd, digBInd, iRow);
                if isnan(prev)
                    prev = 0;
                end
                BDigs(digAInd, digBInd, iRow) = prev + 1;
            end
        end
        
        B = sum(B, 4, 'omitnan');

        B(isnan(B)) = 0;

        BSummed_dim3 = sum(B,3, 'omitnan'); % sum the counts across the angles. the 4x4 matrix will be total counts for each dig pair

        BMean_dim3 = mean(B, 3, 'omitnan'); % this is the average count
        B_normed_dim3 = B ./ BSummed_dim3; % normed by total counts for each dig pair
        B_normed_dim3(~isfinite(B_normed_dim3)) = 0;

        B_normed_dim12 = B ./ sum(B, [1,2], 'omitnan'); % normed by the all counts for a given angle

        R(k).animalName = animalName;
        R(k).sessionName = sessionName;
        R(k).B = B;
        R(k).BDigs = mean(BDigs, 3, 'omitnan');
        R(k).BSummed_dim3 = BSummed_dim3;
        R(k).BMean_dim3 = BMean_dim3;
        R(k).B_normed_dim3 = B_normed_dim3;
        R(k).B_normed_dim12 = B_normed_dim12;
        R(k).rotations = rotations;
        R(k).digTypes = digTypes;
        k = k + 1;
    end % iSession
    
    [SessionsToGroups] = load_sessions_to_groups_table(dataset.analysisSettings);
    R = struct2table(R);
    eliminateUnmatched = true;
    R = helper_add_group_labels_to_table(R, SessionsToGroups, eliminateUnmatched);
end % function



%% FUNCTIONS
function [NormalPerCell90] = recompute_bfo90_percell(dataset)
    isCalciumData = dataset.analysisSettings.IS_CALCIUM_DATA;    
    MapsSquareData = dataset.analysisInput.MapsSquareData;
    analysisSettings = dataset.analysisSettings;
        
    %% Compute the PerCell table used for the NORMAL BFO 90
    % There will be one row per unique cell
    [NormalPerCell90, success] = compute_percell_and_behaviour_for_table(MapsSquareData, analysisSettings.NORMAL_BFO90_ROTATIONS_DEG, analysisSettings.NORMAL_BFO90_CONTEXTS_TO_REFLECT_MAPS, analysisSettings.NORMAL_BFO90_MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    if ~success
        fprintf('Error while processing PerCell correlations\n');
    end

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
