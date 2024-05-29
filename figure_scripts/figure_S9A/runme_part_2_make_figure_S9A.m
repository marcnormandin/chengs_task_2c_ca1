% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This code makes Figure S9A (but used Excel + R for paper).
close all
clear all
clc

% Paths
run('../load_figure_config.m')
INPUT_FOLDER = './local_output';

OUTPUT_FOLDER = CHENGS_TASK_2C_FIGURES_CONFIG.FIGURE_OUTPUT_FOLDER;
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
end

calData = readtable(fullfile(INPUT_FOLDER, 'information_content_calcium_data.xlsx'));


%%
clc
close all

R = calData;
x = prepare_x(R);

fprintf('Plotting\n');
groupLabels = {'Day 1', 'Day 2', 'Day 3'};

hFig = figure('position', get(0,'screensize'));
bpg = boxplotGroup(x, 'primaryLabels', {'FI Cells', 'FS Cells'}, 'secondaryLabels', groupLabels);
legend('location', 'north')
ylabel('Spatial Information Content (bits/sec)')
set(gca, 'fontweight', 'bold', 'fontsize', 18);
a = axis;
axis([a(1), a(2), 0, 100]);

mulana_savefig(hFig, OUTPUT_FOLDER, 'figure_S9A', {'png', 'svg'});

%% Save data to excel as wanted by natcomms
stabilityTypes = {'FI', 'FS'};
nGroups = length(groupLabels);
k = 1;
columnNames = {'A1', 'B1', 'C1', 'D1', 'E1', 'F1'};
for iStability = 1:length(stabilityTypes)
    for iGroup = 1:nGroups
        columnTitle = sprintf('calcium_information_content_%s_%s', stabilityTypes{iStability}, strrep(groupLabels{iGroup}, ' ', '_'));
        X = array2table(x{iStability}(:,iGroup), 'VariableNames', {columnTitle});
        writetable(X, fullfile(OUTPUT_FOLDER, "natcomms_excel_figure_S9A.xlsx"), 'Sheet', 'figure_S9A', 'Range', columnNames{k})
        k = k + 1;
    end
end % iStability

function [x] = prepare_x(R)
    stabilities = [true, false];
    nStabilities = length(stabilities);
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    nGroups = length(groupLabels);

    x = cell(1, nStabilities);
    disp(size(x))
    fprintf('Organizing\n');
    for iStability = 1:nStabilities
        stability = stabilities(iStability);

        xDataCell = cell(1,nGroups);
        for iGroup = 1:nGroups
            groupLabel = groupLabels{iGroup};
            xDataCell{iGroup} = R.ics_avg(ismember(R.groupLabel, groupLabel) & R.isStable == stability);
        
            xDataCell{iGroup} = rmoutliers(xDataCell{iGroup}, 'percentiles', [0, 95]);
        end % iGroup

        % nLargest
        nLargest = -1;
        for iGroup = 1:nGroups
            nMembers = length(xDataCell{iGroup});
            if nMembers > nLargest
                nLargest = nMembers;
            end
        end
        xData = nan(nLargest, nGroups);
        for iGroup = 1:nGroups
            nMembers = length(xDataCell{iGroup});
            xData(1:nMembers, iGroup) = xDataCell{iGroup};
        end
        x{iStability} = xData;
    end % iGroup
end % function

