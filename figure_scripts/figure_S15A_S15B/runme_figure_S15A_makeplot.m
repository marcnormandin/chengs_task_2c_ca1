% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This function reads in the excel files and then plots them since they
% were previously plotted in Excel. This is so that a reviewer doesn't have
% to plot them manually. Makes Figure S15A (mfr).

close all
clear all
clc

% Paths
run('../load_figure_config.m')
OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;


% This folder should have subfolders for each occupancy threshold.
MAIN_FOLDER = fullfile('.\local_output','mfrs');
FRACTIONS = [0.5, 0.6, 0.7];

% To make the grouped boxplot I need to modify the table.
fnList = cell(3,2);
for i = 1:length(FRACTIONS)
    fraction = FRACTIONS(i);

    fnList{i,1} = fraction;
    fnList{i,2} = fullfile(MAIN_FOLDER, sprintf('min_fraction_occupied_%0.1f', fraction), sprintf('tetrodes_rate_matrix_within_across_per_cell_mfrs_min_fraction_occupied_%0.1f.xlsx', fraction));
end

hFig = figure('position', get(0, 'screensize'));
for i = 1:size(fnList,1)
    fraction = fnList{i,1};
    fn = fnList{i,2};
    x = prepare_x(fn);

    subplot(1,size(fnList,1), i);
    boxplotGroup(x,'primaryLabels',{'within', 'across'}, 'secondaryLabels', {'Day 1', 'Day 2', 'Day 3'})
    a = axis;
    axis([a(1), a(2), 0, 0.4])
    subtitle(sprintf('Fraction %0.1f', fraction))
    %set(gca, 'fontweight', 'bold', 'fontsize', 18)
end

mulana_savefig(hFig, OUTPUT_FOLDER, "figure_S15A", {'png', 'svg'});


%% Added to export to excel for natcomms
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroups = length(groupLabels);
excelColumns = string(('a':'z').').';
k = 1;
for i = 1:size(fnList,1)
    fraction = fnList{i,1};
    fn = fnList{i,2};
    x = prepare_x(fn);

    comparisonTypes = {'within', 'across'};
    for iComparison = 1:length(comparisonTypes)
        comparisonType = comparisonTypes{iComparison};
        for iGroup = 1:nGroups
            dayName = groupLabels{iGroup};

            columnName = sprintf('mfrs_fraction_%0.1f_%s_%s', fraction, comparisonType, strrep(groupLabels{iGroup}, ' ', '_'));
            xSave = x{iComparison}(:,iGroup);
            
            writetable(array2table(xSave, 'VariableNames', {columnName}), fullfile(OUTPUT_FOLDER, "natcomms_excel_figure_S15A.xlsx"), 'Sheet', 'figure_S15A', 'Range', [excelColumns(k) '1']);
            
            k = k + 1;
        end
    end % iComparison
end



%% Functions
function [x] = prepare_x(fn)
    T = readtable(fn);
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

