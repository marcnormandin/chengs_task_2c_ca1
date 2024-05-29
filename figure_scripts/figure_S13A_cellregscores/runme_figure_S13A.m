% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This code makes Figure S13A (cellreg scores)

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
% LOCAL_OUTPUT_FOLDER = './local_output';
% if ~exist(LOCAL_OUTPUT_FOLDER, 'dir')
%     mkdir(LOCAL_OUTPUT_FOLDER);
% end

calciumAnalysisResultsFolder = fullfile(INPUT_FOLDER, 'calcium_20220516_170618');

%%
calData = load_dataset(calciumAnalysisResultsFolder);
%%
BestAligned = calData.BestAligned;
tmp = load(fullfile(INPUT_FOLDER, "calcium_CellRegTable.mat"));
CellRegTable = tmp.CellRegTable;
clear tmp

% Add the scores as a column of the BestAligned table
nRows = size(BestAligned,1);
cellScores = nan(nRows,1);
cellRegNSessionsPresent = nan(nRows,1);
for iRow = 1:nRows
    cellName = BestAligned.cellName{iRow};
    indMatch = find(ismember(CellRegTable.cellName, cellName));
    if length(indMatch) ~= 1
        error('Something is inconsistent.\n');
    end
    cellScores(iRow) = CellRegTable.cellScore(indMatch);
    cellRegNSessionsPresent(iRow) = CellRegTable.sharedCount(indMatch);
end
BestAligned.cellRegScore = cellScores;
BestAligned.cellRegNSessionsPresent = cellRegNSessionsPresent;
fprintf('Done adding the scores.\n');

%% Deciles
BA = BestAligned;
BA(BA.cellRegNSessionsPresent < 2,:) = []; % Cells must be present in 2+ sessions.
% nBins = 10; % deciles
% edges = linspace(0,1,nBins+1);
% p = prctile(BA.bestCorrelation, edges(2:end-1)*100);
% pe = [0, p, 1];
p = [0.3, 0.42, 0.51, 0.58, 0.64, 0.69, 0.75, 0.80, 0.86];
BA.binId = discretize(BA.bestCorrelation, [0, p, 1]);
BA.binId(BA.bestCorrelation < 0) = 1;
% By def nan is FS
BA.binId(isnan(BA.bestCorrelation)) = 1;

close all
hFig = figure();
boxplot(BA.cellRegScore, BA.binId)
xlabel('Decile')
ylabel('Cell Registration Score')
title('Cell Registration by Decile')
grid on
set(gca, 'fontweight', 'bold', 'fontsize', 18);

mulana_savefig(hFig, OUTPUT_FOLDER, "figure_S13A", {'png', 'svg'});

%% Added to export data for natcomms
excelRanges = {'A1', 'B1', 'C1', 'D1', 'E1', 'F1','G1','H1','I1','J1'};
for binId = 1:10
    columnName = sprintf('cellregscore_decile_%d', binId);
    x = BA.cellRegScore(BA.binId == binId);
    X = array2table(x, 'VariableNames', {columnName});
    writetable(X, fullfile(OUTPUT_FOLDER, "natcomms_excel_figure_S13A.xlsx"), 'Sheet', 'figure_S13A', 'Range', excelRanges{binId});
end
