% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This reproduces Figure 6B (which was made in Excel).

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
LOCAL_OUTPUT_FOLDER = './local_output'; % for Figure S5B
if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
    mkdir(LOCAL_OUTPUT_FOLDER);
end

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');
calResultsData = load(fullfile(calciumAnalysisResultsFolder, 'analysis_results.mat'));
calInputData = load(fullfile(calciumAnalysisResultsFolder, 'analysis_input.mat'));

%% We need the BestAligned table and the CellRegTable
BestAligned = calResultsData.analysisResults.BestAligned;
CellRegTable = calInputData.analysisInput.CellRegTable;
SessionsToGroups = load_calcium_sessions_to_groups_table();

BestAligned = helper_add_group_labels_to_table(BestAligned, SessionsToGroups, true);

%% Add the registered cell names to the BestAligned table
registeredCellNames = cell(size(BestAligned,1),1);
for iRow = 1:size(BestAligned,1)
    sessionCellName = BestAligned.cellName{iRow};
    % Search the CellRegTable for the same now and then get the registered
    % name
    indMatch = find(ismember(CellRegTable.cellName, sessionCellName));
    if length(indMatch) ~= 1
        error('Something is wrong.');
    end
    registeredCellName = CellRegTable.registeredCellName{indMatch};

    registeredCellNames{iRow} = registeredCellName;
end
BestAligned.registeredCellName = registeredCellNames;

%%
BA = BestAligned;
groupLabels = {'Day 1', 'Day 2', 'Day 3'};

%pairs = nchoosek(groupLabels, 2);
pairs = cell(3,2);
pairs{1,1} = 'Day 1';
pairs{1,2} = 'Day 2';

pairs{2,1} = 'Day 2';
pairs{2,2} = 'Day 3';

pairs{3,1} = 'Day 1';
pairs{3,2} = 'Day 3';

nPairs = size(pairs,1);
R = [];
rk= 1;
for iPair = 1:nPairs
    groupLabelA = pairs(iPair,1);
    groupLabelB = pairs(iPair,2);

    BAA = BA(ismember(BA.groupLabel, groupLabelA),:);
    regCellsA = unique(BAA.registeredCellName);

    BAB = BA(ismember(BA.groupLabel, groupLabelB),:);
    regCellsB = unique(BAB.registeredCellName);

    % Find the ones common to both days
    regCommon = intersect(regCellsA, regCellsB);

    THRESHOLD = 0.3;

    correlationPairs = nan(length(regCommon),2); % two for the two days
    for iCommon = 1:length(regCommon)
        regCellName = regCommon{iCommon};

        ind1 = find(ismember(BAA.registeredCellName, regCellName));
        ind2 = find(ismember(BAB.registeredCellName, regCellName));

        if isempty(ind1) || isempty(ind2)
            error('Something is wrong');
        end

        correlationPairs(iCommon, 1) = BAA.bestCorrelation(ind1);
        correlationPairs(iCommon, 2) = BAB.bestCorrelation(ind2);
    end
    nCellsTotalBeforeNanCriteria = size(correlationPairs,1);

    correlationPairs(any(isnan(correlationPairs),2),:) = [];

    nCellsTotalAfterNanCriteria = size(correlationPairs,1);

    stabilityPairs = correlationPairs > THRESHOLD;

    stabilityTypes = {true, false};
    for iStability1 = 1:length(stabilityTypes)
        stability1 = stabilityTypes{iStability1};
        for iStability2 = 1:length(stabilityTypes)
            stability2 = stabilityTypes{iStability2};

            % Now for the pair of days AND the two stability types, get the
            % counts by check each cell's data
            count = 0;
            for iCell = 1:size(stabilityPairs,1)
                if stabilityPairs(iCell,1) == stability1 && stabilityPairs(iCell,2) == stability2
                    count = count + 1;
                end
            end

            R(rk).iPair = iPair;
            R(rk).groupLabelA = groupLabelA;
            R(rk).groupLabelB = groupLabelB;
            R(rk).nCells = length(regCommon);
            R(rk).stability1 = stability1;
            R(rk).stability2 = stability2;
            R(rk).count = count;
            R(rk).nCellsTotalBeforeNanCriteria = nCellsTotalBeforeNanCriteria;
            R(rk).nCellsTotalAfterNanCriteria = nCellsTotalAfterNanCriteria;
            rk = rk + 1;
        end
    end
end
R = struct2table(R);

%% Now go through to calculate the percents
clc
P = [];
pk = 1;
stabilities = [true, false];
for iPair = 1:nPairs
    groupLabelA = pairs(iPair,1);
    groupLabelB = pairs(iPair,2);

    for iStability = 1:length(stabilities)
        stability = stabilities(iStability);

        % Get the rows that match the day pair and the day A stability type
        M = R(ismember(R.groupLabelA, groupLabelA) & ismember(R.groupLabelB, groupLabelB) & R.stability1 == stability,:);

        nTotal = sum(M.count);
        nSame = sum(M.count(M.stability2 == stability));
        nDifferent = nTotal - nSame;
        if nDifferent ~= sum(M.count(M.stability2 ~= stability))
            error('Bad math');
        end
        P(pk).iPair = iPair;
        P(pk).groupLabelA = groupLabelA;
        P(pk).groupLabelB = groupLabelB;
        P(pk).stability = stability;
        P(pk).nTotal = nTotal;
        P(pk).nSame = nSame;
        P(pk).nDifferent = nDifferent;
        P(pk).percentSame = nSame / nTotal * 100;
        P(pk).percentDifferent = nDifferent / nTotal * 100;
        pk = pk + 1;
    end % iStability
end % iPair
P = struct2table(P);


%% Make the figure
close all
clc

hFig = figure('position', get(0, 'screensize'));
PP = 2;
QQ = 3;
kk = 1;
ax = [];
stabilityNames = {'FI', 'FS'};
for iStability = 1:length(stabilities)
    stability = stabilities(iStability);
    stabilityName = stabilityNames{iStability};
    for iPair = 1:nPairs
        groupLabelA = pairs{iPair,1};
        groupLabelB = pairs{iPair,2};

        ax(kk) = subplot(PP,QQ,kk);
        kk = kk + 1;
        plot_pie(P, groupLabelA, groupLabelB, stability)
        title(sprintf('%s %s -> %s', stabilityName, groupLabelA, groupLabelB));
    end
end

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_6B_piecharts', {'png', 'svg'});

%% Added to export for NatComms
% Use the FI and FS nomenclature
NC = P;
nRows = size(NC,1);
stabilityString = cell(nRows,1);
for iRow = 1:nRows
    if NC.stability(iRow) == true
        stabilityString{iRow} = 'FI';
    else
        stabilityString{iRow} = 'FS';
    end
end
NC.type = stabilityString;
NC(:, ismember(NC.Properties.VariableNames, {'iPair', 'stability'})) = [];
writetable(NC, fullfile(OUTPUT_FOLDER, 'natcomms_excel_figure_6B.xlsx'), 'Sheet', 'figure_6B')



function plot_pie(P, groupLabelA, groupLabelB, stability)
    X = P(ismember(P.groupLabelA, groupLabelA) & ismember(P.groupLabelB, groupLabelB) & P.stability == stability, :);
    x = [X.percentSame(1); X.percentDifferent(1)];
    pie(x)
end
