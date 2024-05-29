% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This makes Figure S3. Roseplots for Calcium and Tetrodes.
% It needs to be run twice to make the two figures. Set DO_CALCIUM to false
% for Figure S3B and DO_CALCIUM to true for Figure S3C.
%
% Note that the figures were post-processed to make them have a solid
% color for the interior.
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


%% Load the calcium results from the main analysis since we want the maps which are saved to a mat file.
DO_CALCIUM = false; % Set this to true to process calcium and false for tetrodes.

if DO_CALCIUM
    analysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
    load(fullfile(analysisResultsFolder, 'analysis_input.mat'));
    load(fullfile(analysisResultsFolder, 'analysis_results.mat'));
else
    analysisResultsFolder = fullfile(INPUT_FOLDER, 'tetrodes_20220316_100902');
    load(fullfile(analysisResultsFolder, 'analysis_input.mat'));
    load(fullfile(analysisResultsFolder, 'analysis_results.mat'));
end


%%
if analysisSettings.IS_CALCIUM_DATA
    datasetName = 'calcium';
else
    datasetName = 'tetrodes';
end

%%
close all
clc

REFLECT_SECOND_CONTEXT = false;
ROTATED_MAPS_BASED_ON_BEST_ALIGNED_SEQUENCE = false;

MapsData = analysisInput.MapsData;
BestAligned = analysisResults.BestAligned;

% We want to operate over cells
[T] = get_sessions_and_cellnames_from_table(MapsData);

% Results
R = [];
k = 0;
for iRow = 1:size(T,1)
    animalName = T.animalName{iRow};
    sessionName = T.sessionName{iRow};
    cellName = T.cellName{iRow};
    animalSessionCellName = T.animalSessionCellName{iRow};
    
    fprintf('Processing %d of %d: %s %s %s\n', iRow, size(T,1), animalName, sessionName, cellName);
    
    % Get all of the data associated with the cell (al trials)
    cellData = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);

    baCellData = BestAligned(ismember(BestAligned.animalSessionCellName, cellData.animalSessionCellName),:);
    if isempty(baCellData)
        error('No best aligned data found for %s', cellData.animalSessionCellName);
    end
    if ~all(baCellData.cellDataInds{1} == cellData.cellDataInds)
        error('these should be 1 to 1');
    end


    % Get the best aligned row
    if ROTATED_MAPS_BASED_ON_BEST_ALIGNED_SEQUENCE
        rotatedCellDataInds = [];
    
        c1inds = find(cellData.contextIds == 1);
        c1CellDataInds = cellData.cellDataInds(c1inds);
        r1 = find(baCellData.rotationSequence1{1} == 1);
        rotatedCellDataInds = cat(1, rotatedCellDataInds(:), c1CellDataInds(r1));
    
    
        c2inds = find(cellData.contextIds == 2);
        c2CellDataInds = cellData.cellDataInds(c2inds);
        r2 = find(baCellData.rotationSequence2{1} == 1);
        rotatedCellDataInds = cat(1, rotatedCellDataInds(:), c2CellDataInds(r2));

        disp(r1)
        disp(r2)
        disp(rotatedCellDataInds)
    
        for iMap = 1:cellData.numEntries
            if any(ismember(rotatedCellDataInds, cellData.cellDataInds(iMap)))
                % If it is one of the maps that need to be rotated, then rotate
                % it
                m = cellData.maps(:,:,iMap);
                m = rot90(m,2);
                cellData.maps(:,:,iMap) = m; % store it back
                fprintf('Rotating map\n');
            end
        end
    end % if
    
    % Remove any empty maps that we don't want to use
    badinds = find(all(isnan(cellData.maps), [1,2]));
    if ~isempty(badinds)
        fprintf('Found %d empty maps.\n', length(badinds));

        % Eliminate 
        cellData.trialIds(badInds) = [];
        cellData.contextIds(badInds) = [];
        cellData.digs{badInds} = [];
        cellData.cellDataInds(badInds) = [];
        cellData.maps(:,:,badinds) = [];
    end

    % Surviving trials
    nTrials = length(cellData.trialIds);
    
    % A
    for i = 1:nTrials
        % B
        for j = i+1:nTrials
            if i == j
                continue; % skip map with itself
            end
            mapA = cellData.maps(:,:,i);
            mapB = cellData.maps(:,:,j);

            contextA = cellData.contextIds(i);
            contextB = cellData.contextIds(j);

            if contextA ~= contextB
                % across comparison
                if contextA ~= 1
                    tmp = mapA;
                    mapA = mapB;
                    mapB = tmp;
                    
                end

                if REFLECT_SECOND_CONTEXT
                    mapB = fliplr(mapB);
                end
            end

            if all(mapB==0, 'all') || all(mapA==0, 'all')
                continue;
            end
            
            % Compute the center out angles and the locations
            [aAngle, ax, ay] = ml_alg_placemap_center_out_angle(mapA, []);
            [bAngle, bx, by] = ml_alg_placemap_center_out_angle(mapB, []);
            % Compute the angle difference 
            [angleDiff] =  ml_alg_center_out_difference(aAngle, bAngle);

            % Store
            k = k + 1;
            R(k).animalName = animalName;
            R(k).sessionName = sessionName;
            R(k).cellName = cellName;
            R(k).animalSessionCellName = animalSessionCellName;
            R(k).cellId = cellData.cellId;
            R(k).bestCorrelation = baCellData.bestCorrelation;
            
            % A
            R(k).trialId_a = cellData.trialIds(i);
            R(k).contextId_a = cellData.contextIds(i);
            R(k).cellDataInd_a = cellData.cellDataInds(i);
            R(k).dig_a = cellData.digs{i};
            R(k).centerOutAngle_a = aAngle;
            %R(k).dirv = dirvA;

            R(k).centerOutLocation_a = [ax, ay];

            % B
            R(k).trialId_b = cellData.trialIds(j);
            R(k).contextId_b = cellData.contextIds(j);
            R(k).cellDataInd_b = cellData.cellDataInds(j);
            R(k).dig_b = cellData.digs{j};
            R(k).centerOutAngle_b = bAngle;
            %R(k).dirv = dirvB;

            R(k).centerOutLocation_b = [bx, by];

            if cellData.contextIds(i) == cellData.contextIds(j)
                R(k).comparisonClassification = 'within';
            else
                R(k).comparisonClassification = 'across';
            end
            R(k).centerOutAngleDifference = angleDiff;
        end
    end
end

R = struct2table(R);

% Add the group/days to the table
R = helper_add_group_labels_to_table(R, load_sessions_to_groups_table(analysisSettings), true);

%% Plot only per day (nothing else matters, stability, within/across)
close all
clc

RTICK_VALUES = 0:0.01:0.09;
LL = max(RTICK_VALUES);

groupLabels = {'Day 1', 'Day 2', 'Day 3'}; % unique(R.groupLabel);
nGroups = length(groupLabels);
nRows = 1;
nCols = nGroups;

hFig = figure('position', get(0, 'screensize'));
for iGroup = 1:nGroups
    groupLabel = groupLabels{iGroup};

    S = R(ismember(R.groupLabel, groupLabel),:);
    r = [];
    r.coa_deg = S.centerOutAngleDifference(:);

    if ~isempty(r.coa_deg)
        % Plot
        subplot(nRows,nCols,iGroup);

        fprintf('Processing %s:\n', groupLabel);

        % mean about 0
        fprintf('About 0 degrees\n');
        [mu0, ll0, ul0] = compute_circ_mean_stats_0(r.coa_deg);
        polarplot([deg2rad(mu0), deg2rad(mu0)], [0, LL], 'r-', 'linewidth', 2);
        hold on
        polarplot([deg2rad(ll0), deg2rad(ll0)], [0, LL], 'r--', 'linewidth', 2);
        polarplot([deg2rad(ul0), deg2rad(ul0)], [0, LL], 'r--', 'linewidth', 2);
        set(gca, 'rtick', RTICK_VALUES)
        set(gca, 'rlim', [0, max(RTICK_VALUES)]);

        % mean about 180
        fprintf('About 180 degrees\n');
        [mu180, ll180, ul180] = compute_circ_mean_stats_180(r.coa_deg);
        polarplot([deg2rad(mu180), deg2rad(mu180)], [0, LL], 'm-', 'linewidth', 2);
        polarplot([deg2rad(ll180), deg2rad(ll180)], [0, LL], 'm--', 'linewidth', 2);
        polarplot([deg2rad(ul180), deg2rad(ul180)], [0, LL], 'm--', 'linewidth', 2);
        set(gca, 'rtick', RTICK_VALUES)
        set(gca, 'rlim', [0, max(RTICK_VALUES)]);
        fprintf('\n\n');

        % Main roseplot
        polarplot_angle_differences(r.coa_deg, 0:12:360);

        title(sprintf('%s\nn =%d\n%0.2f [%0.2f, %0.2f]\n%0.2f [%0.2f, %0.2f]', groupLabel, length(r.coa_deg), mu0, ll0, ul0, mu180, ll180, ul180), 'interpreter', 'none');
    end
end % iGroup

ml_savefig(hFig, OUTPUT_FOLDER, sprintf('figure_S3_%s_roseplots_per_day', datasetName), {'png', 'svg', 'fig'});

%% Added to save to excel for natcomms
LOCAL_OUTPUT_FOLDER = fullfile(pwd, 'local_output');
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

%For Tetrodes ABC, for Calcium DEF
rangeLetters = 'ABCDEF';
groupLabels = {'Day 1', 'Day 2', 'Day 3'}; % unique(R.groupLabel);
nGroups = length(groupLabels);
for iGroup = 1:nGroups
    groupLabel = groupLabels{iGroup};
    S = R(ismember(R.groupLabel, groupLabel),:);

    columnName = sprintf('center_out_angle_difference_%s_%s', datasetName, strrep(groupLabel, ' ', '_'));
    X = S.centerOutAngleDifference(:);
    XT = array2table(X,'VariableNames',{columnName});
    nRows = size(XT,1)+1; % +1 for the column name
    writetable(XT, fullfile(LOCAL_OUTPUT_FOLDER, 'natcomms_excel_figure_S3.xlsx'), 'Sheet', 'figure_S3', 'Range', [rangeLetters(DO_CALCIUM*nGroups+iGroup) '1'])

end % iGroup


function polarplot_angle_differences(angles, binEdges)    
    hc = histcounts(angles, binEdges);
    p = hc ./ sum(hc);

    binCenters = binEdges(1:end-1) + diff(binEdges)/2;
    theta = deg2rad(binCenters); 

    % Connect the final to the first point
    theta = [theta, theta(1)];
    rho = [p, p(1)];

    polarplot(theta, rho, 'b-', 'linewidth', 4);
end


function [a, p] = get_pmf(angles, binEdges)
    a = binEdges(1:end-1) + diff(binEdges)/2;
    hc = histcounts(angles, binEdges);
    p = hc ./ sum(hc);
end

function [mu, ll, ul] = compute_circ_mean_stats_0(x)
    x0 = x(x >= 270 | x <= 90); % peak around zero
    [mu, ll, ul] = compute_circ_mean_stats(x0);
end

function [mu, ll, ul] = compute_circ_mean_stats_180(x)
    x180 = x(x >= 90 & x <= 270); % peak around 180
    [mu, ll, ul] = compute_circ_mean_stats(x180);
end

function [mu, ll, ul] = compute_circ_mean_stats(y)
    [mu, ul, ll] = circ_mean(deg2rad(y));
    
    mu = rad2deg(mu);
    mu(mu <0) = mu + 360;
    
    ul = rad2deg(ul);
    ul(ul < 0) = ul + 360;
    
    ll = rad2deg(ll);
    ll(ll < 0) = ll + 360;
    
    fprintf('The mean is %0.1f degrees (%0.2f, %0.2f)\n', mu, ll, ul);
end % function
