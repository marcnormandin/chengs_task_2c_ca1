% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This part makes the Excel file which is then plotted in part 3.
close all
clear all
clc

SPEED_THRESHOLD_CM_PER_S = 2;
SPIKE_THRESHOLD = 0;

% These ones include m4 which is the PFR
LOCAL_INPUT_FOLDER = './local_output';
OUTPUT_FOLDER = './local_output';

groupLabels = {'Day 1', 'Day 2', 'Day 3'};

T = [];
files = dir(fullfile(LOCAL_INPUT_FOLDER, '*_rm.mat'));
for iFile = 1:length(files)
    fn = fullfile(files(iFile).folder, files(iFile).name);
    tmp = load(fn);
    if isempty(T)
        T = tmp.RMA;
    else
        T = [T; tmp.RMA];
    end
end

RM_ROWS = tmp.RM_ROWS;
NUM_TRIALS = tmp.NUM_TRIALS;

T.groupName = repmat('feature_rich', size(T,1),1);
T = helper_add_group_labels_to_table(T, load_calcium_sessions_to_groups_table(), true);

RateMatricesWithAcrossPerSession = compute_per_animal_session_rate_differences_within_across(T);
RateMatricesWithAcrossPerSession(~ismember(RateMatricesWithAcrossPerSession.groupLabel, groupLabels),:) = [];

writetable(RateMatricesWithAcrossPerSession, fullfile(OUTPUT_FOLDER, 'calcium_rate_matrices_within_across_per_animal_session.xlsx'))


%%
% I copied this from the main analysis so it is consistent with what
% we used for the tetrodes.
function [R] = compute_per_animal_session_rate_differences_within_across(T)
    fieldName = 'm3';
    numRows = size(T,1);
    rateMatrixContextIds = [zeros(1,6), ones(1,6)];
    R = [];
    k = 1;
    for iRow = 1:numRows        
        groupId = T.groupId(iRow);
        groupLabel = T.groupLabel{iRow};
        rateMatrixMean = T.(fieldName){iRow}; 

        withinData = [];
        acrossData = [];
        for i = 1:length(rateMatrixContextIds)
            for j = 1+1:length(rateMatrixContextIds)
                x = rateMatrixMean(i,j);
                if rateMatrixContextIds(i) ~= rateMatrixContextIds(j)
                    acrossData = cat(1, acrossData, x);
                else
                    withinData = cat(1, withinData, x);
                end
            end
        end

        withinMean = nanmean(withinData, 'all');
        withinError = nanstd(withinData, 1, 'all') ./ sqrt(length(withinData));

        acrossMean = nanmean(acrossData, 'all');
        acrossError = nanstd(acrossData, 1, 'all') ./ sqrt(length(acrossData));

        R(k).animalName = T.animalName{iRow};
        R(k).sessionName = T.sessionName{iRow};
        R(k).groupId = groupId;
        R(k).groupLabel = groupLabel;
        R(k).withinMean = withinMean;
        R(k).withinError = withinError;
        R(k).acrossMean = acrossMean;
        R(k).acrossError = acrossError;

        k = k + 1;
    end % iRow
    R = struct2table(R);
end % function
