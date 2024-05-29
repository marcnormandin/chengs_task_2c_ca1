% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Figure S16A and Figure S16B (Absolute Rate Difference for Tetrodes separated by cell type)
% This uses the min occupancy of 50%.

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

fn = '..\figure_S15A_S15B\local_output\mfrs\min_fraction_occupied_0.5\tetrodes_rate_matrix_within_across_per_cell_mfrs_min_fraction_occupied_0.5.xlsx';

T = readtable(fn);


%% Box Plot
hFig = figure('position', get(0, 'screensize'));
ax = [];

% A) Stable / FI Cells
ax(1) = subplot(1, 2, 1);
x = prepare_x(T, 'stable');
boxplotGroup(x,'primaryLabels',{'within', 'across'}, 'secondaryLabels', {'Day 1', 'Day 2', 'Day 3'})
a = axis;
axis([a(1), a(2), 0, 0.6])
ylabel('Absolute Rate Difference (Hz)');
title('FI Cells')

% B) Unstable / FS Cells
ax(1) = subplot(1, 2, 2);
x = prepare_x(T, 'unstable');
boxplotGroup(x,'primaryLabels',{'within', 'across'}, 'secondaryLabels', {'Day 1', 'Day 2', 'Day 3'})
a = axis;
axis([a(1), a(2), 0, 0.6])
ylabel('Absolute Rate Difference (Hz)');
title('FS Cells')

linkaxes(ax, 'xy')

mulana_savefig(hFig, CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER, 'figure_S16A_S16B', {'png', 'svg'});




%% Added to export to excel for natcomms
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroups = length(groupLabels);
excelColumns = string(('a':'z').').';
k = 1;
stabilityTypes = {'stable', 'unstable'};
cellTypes = {'FI', 'FS'};
for iStabilityType = 1:length(stabilityTypes)
    stabilityType = stabilityTypes{iStabilityType};
    cellType = cellTypes{iStabilityType};

    x = prepare_x(T, stabilityType);
    
    comparisonTypes = {'within', 'across'};
    for iComparison = 1:length(comparisonTypes)
        comparisonType = comparisonTypes{iComparison};
        for iGroup = 1:nGroups
            dayName = groupLabels{iGroup};
    
            columnName = sprintf('absolute_rate_difference_%s_%s_%s', cellType, comparisonType, strrep(groupLabels{iGroup}, ' ', '_'));
            xSave = x{iComparison}(:,iGroup);
            
            writetable(array2table(xSave, 'VariableNames', {columnName}), fullfile(OUTPUT_FOLDER, "natcomms_excel_figure_S16A_S16B.xlsx"), 'Sheet', 'figure_S16A_S16B', 'Range', [excelColumns(k) '1']);
            
            k = k + 1;
        end
    end % iComparison
end % iStabilityType





%% Functions

function [x] = prepare_x(T, stabilityType)
    % stability type is stable or unstable
    T = T(ismember(T.stability, stabilityType),:); % filter for the type

    badInds = union(find(isnan(T.within_mean)), find(isnan(T.across_mean)));
    T(badInds,:) = [];
    
    xw1 = T.within_mean(T.dayId == 1);
    nw1 = length(xw1);
    xw2 = T.within_mean(T.dayId == 2);
    nw2 = length(xw2);
    xw3 = T.within_mean(T.dayId == 3);
    nw3 = length(xw3);
    nw = max([nw1, nw2, nw3]);
    xw = nan(nw, 3);
    xw(1:nw1,1) = xw1;
    xw(1:nw2,2) = xw2;
    xw(1:nw3,3) = xw3;
    
    xa1 = T.across_mean(T.dayId == 1);
    na1 = length(xa1);
    xa2 = T.across_mean(T.dayId == 2);
    na2 = length(xa2);
    xa3 = T.across_mean(T.dayId == 3);
    na3 = length(xa3);
    na = max([na1, na2, na3]);
    xa = nan(na, 3);
    xa(1:na1,1) = xa1;
    xa(1:na2,2) = xa2;
    xa(1:na3,3) = xa3;
    x = {xw, xa};
end % function


function [y] = make_y(T, stabilityType)
    y = nan(3,2);
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    columns = {'within_mean', 'across_mean'};
    for iGroup = 1:length(groupLabels)
        groupLabel = groupLabels{iGroup};
        for iColumn = 1:length(columns)
            columnName = columns{iColumn};
    
            y(iGroup, iColumn) = mean(T.(columnName)(ismember(T.dayLabel, groupLabel) & ismember(T.stability, stabilityType)), 'omitnan');
        end % iColumn
    end % iGroup
end % function

