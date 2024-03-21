% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
%
function compute_popvectors_across_contexts(analysisSettings, CellRegTable, BestAligned, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, STABILITY_CLASSIFICATION_GROUP_LABEL)

%%
% Add the group labels (day names) to the data
SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
eliminateUnmatched = true;
BestAligned = helper_add_group_labels_to_table(BestAligned, SessionsToGroups, eliminateUnmatched);
fprintf('Added the day/group label information to the BestAligned results.\n');

% TEMP
% BestAligned(BestAligned.metricValue1 <= 0.5 | BestAligned.metricValue2 <= 0.5,:) = [];
% warning('using a temp, not final!!!');

groupLabels = {'Day 1', 'Day 2', 'Day 3'};

% Classify the stability of the cells based on the given day.
[registeredCellNamesStable, registeredCellNamesUnstable] =  get_stable_unstable_registered_cell_names_classified_1day(analysisSettings, CellRegTable, BestAligned, STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);
fprintf('There are %d stable cells, and %d unstable cells, based on %s with a minimum correlation threshold of %0.1f\n', length(registeredCellNamesStable), length(registeredCellNamesUnstable), STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);

% Add the associated registered cell names to the best aligned table
BestAligned = bestaligned_add_registered_cell_names(BestAligned, CellRegTable);

% We only care about these "Days".
BestAligned(~ismember(BestAligned.groupLabel, groupLabels),:) = [];

% Note, that we might want to update the 'isStable' column of the
% BestAligned table, but the code works without it because we select by
% stable and unstable registered cell names and that used the input
% correlation from the main function input, not that necessarily set by the
% main analysis settings.
BestAlignedStable = BestAligned(ismember(BestAligned.registeredCellName, registeredCellNamesStable),:);
BestAlignedUnstable = BestAligned(ismember(BestAligned.registeredCellName, registeredCellNamesUnstable), :);

% Go through each day and compute the cell type counts
for iGroup = 1:length(groupLabels)
    groupLabel = groupLabels{iGroup};
    T1 = BestAlignedStable(ismember(BestAlignedStable.groupLabel, groupLabel),:);
    T2 = BestAlignedUnstable(ismember(BestAlignedUnstable.groupLabel, groupLabel),:);
    numStable = size(T1,1);
    numUnstable = size(T2,1);
    numTotal = numStable + numUnstable;
    fprintf('Stable cells make up %0.2f of cells on %s\n', numStable/numTotal*100, groupLabel);
end



% Compute stuff with bestaligned maps across days

% Compute the within and across for all combinations of days for both types
% of cells.
[FStable] = compute_F(BestAlignedStable);
[FUnstable] = compute_F(BestAlignedUnstable);

% Isabel said to classify stable and unstable this way.
% Stable cells are stable on Day 1 AND stable on Day 3.
% Unstable cells are unstable on Day 1 AND unstable on Day 3.
sessionNameA = 'Day 1';
sessionNameB = 'Day 3';
FS = FStable(ismember(FStable.sessionNameA, sessionNameA) & ismember(FStable.sessionNameB, sessionNameB),:);
FU = FUnstable(ismember(FUnstable.sessionNameA, sessionNameA) & ismember(FUnstable.sessionNameB, sessionNameB),:);

%

% Store the data as a mat file for Celia to get the curve data for final
% figures
popvectorResults = [];
popvectorResults.type = 'across';
popvectorResults.sessionNameA = sessionNameA;
popvectorResults.sessionNameB = sessionNameB;
popvectorResults.BESTALIGNED_STABILITY_THRESHOLD_CRITERIA = BESTALIGNED_STABILITY_THRESHOLD_CRITERIA;
popvectorResults.STABILITY_CLASSIFICATION_GROUP_LABEL = STABILITY_CLASSIFICATION_GROUP_LABEL;

hFig = figure('position', get(0, 'screensize'));

hold on

% Stable (same day)
LW = 4;
FSz = 24;

[popvectorResults.stable.dp_across_day_A.uz, popvectorResults.stable.dp_across_day_A.cz] = ml_alg_cumdist(FS.('dp_across_day_A'){1});
plot(popvectorResults.stable.dp_across_day_A.uz, popvectorResults.stable.dp_across_day_A.cz, 'b-', 'linewidth', LW)
[popvectorResults.stable.dp_across_day_B.uz, popvectorResults.stable.dp_across_day_B.cz] = ml_alg_cumdist(FS.('dp_across_day_B'){1});
plot(popvectorResults.stable.dp_across_day_B.uz, popvectorResults.stable.dp_across_day_B.cz, 'b--', 'linewidth', LW)

% Unstable (same day)
[popvectorResults.unstable.dp_across_day_A.uz, popvectorResults.unstable.dp_across_day_A.cz] = ml_alg_cumdist(FU.('dp_across_day_A'){1});
plot(popvectorResults.unstable.dp_across_day_A.uz, popvectorResults.unstable.dp_across_day_A.cz, 'r-', 'linewidth', LW)
[popvectorResults.unstable.dp_across_day_B.uz, popvectorResults.unstable.dp_across_day_B.cz] = ml_alg_cumdist(FU.('dp_across_day_B'){1});
plot(popvectorResults.unstable.dp_across_day_B.uz, popvectorResults.unstable.dp_across_day_B.cz, 'r--', 'linewidth', LW)

xticks(0:0.2:1)
yticks(0:0.2:1);
set(gca, 'fontsize', FSz)
axis equal
% % Stable (context)
% [popvectorResults.stable.dp_across_context_1.uz, popvectorResults.stable.dp_across_context_1.cz] = ml_alg_cumdist(FS.('dp_across_context_1'){1});
% plot(popvectorResults.stable.dp_across_context_1.uz, popvectorResults.stable.dp_across_context_1.cz, 'g-', 'linewidth', 2)
% 
% [popvectorResults.stable.dp_across_context_2.uz, popvectorResults.stable.dp_across_context_2.cz] = ml_alg_cumdist(FS.('dp_across_context_2'){1});
% plot(popvectorResults.stable.dp_across_context_2.uz, popvectorResults.stable.dp_across_context_2.cz, 'm-', 'linewidth', 2)
% 
% % Unstable (same day)
% [popvectorResults.unstable.dp_across_context_1.uz, popvectorResults.unstable.dp_across_context_1.cz] = ml_alg_cumdist(FU.('dp_across_context_1'){1});
% plot(popvectorResults.unstable.dp_across_context_1.uz, popvectorResults.unstable.dp_across_context_1.cz, 'g--', 'linewidth', 2)
% 
% [popvectorResults.unstable.dp_across_context_2.uz, popvectorResults.unstable.dp_across_context_2.cz] = ml_alg_cumdist(FU.('dp_across_context_2'){1});
% plot(popvectorResults.unstable.dp_across_context_2.uz, popvectorResults.unstable.dp_across_context_2.cz, 'm--', 'linewidth', 2)





% [ks1,p1] = kstest2(FS.('dp_across_day_A'){1},FS.('dp_across_day_B'){1});
% [ks2,p2] = kstest2(FU.('dp_across_day_A'){1},FU.('dp_across_day_B'){1});

% These compare stable with stable and unstable with unstable
[ks1,p1] = kstest2(FS.('dp_across_day_A'){1},FS.('dp_across_day_B'){1});
[ks2,p2] = kstest2(FU.('dp_across_day_A'){1},FU.('dp_across_day_B'){1});

% These compare stable with unstable
[ks3,p3] = kstest2(FS.('dp_across_day_A'){1},FU.('dp_across_day_A'){1});
[ks4,p4] = kstest2(FS.('dp_across_day_B'){1},FU.('dp_across_day_B'){1});
[ks5,p5] = kstest2(FS.('dp_across_day_A'){1},FU.('dp_across_day_B'){1});
[ks6,p6] = kstest2(FS.('dp_across_day_B'){1},FU.('dp_across_day_A'){1});


legend({...
    sprintf('Stable (on %s) context 1 with context 2 on %s',STABILITY_CLASSIFICATION_GROUP_LABEL, sessionNameA), ...
    sprintf('Stable (on %s) context 1 with context 2 on %s',STABILITY_CLASSIFICATION_GROUP_LABEL, sessionNameB), ...
    sprintf('Unstable (on %s) context 1 with context 2 on %s', STABILITY_CLASSIFICATION_GROUP_LABEL, sessionNameA), ...
    sprintf('Unstable (on %s) context 1 with context 2 on %s',STABILITY_CLASSIFICATION_GROUP_LABEL, sessionNameB)}, ...
    'location', 'southoutside');
grid on
axis square
%title(sprintf('Population Vector Dot Products\n KS-test: p unstable (%s with %s) = %0.2f\nKS-test: p stable (%s with %s) = %0.2f', sessionNameA, sessionNameB, p2, sessionNameA, sessionNameB, p1))
xlabel('Dot Product', 'fontsize', FSz)
ylabel('Cumulative', 'fontsize', FSz)

% save the output to here
outputFolder = fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'popvectors');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% save the figure
mulana_savefig(hFig, outputFolder, sprintf('calcium_popvectors_across_context_classified_%s',STABILITY_CLASSIFICATION_GROUP_LABEL), {'png', 'svg', 'fig'});

% Save the data
outputFilename = fullfile(outputFolder, sprintf('calciumPopvectorResults_across_context_%s.mat',STABILITY_CLASSIFICATION_GROUP_LABEL) );
save(outputFilename, 'popvectorResults', 'analysisSettings', ...
    'sessionNameA', 'sessionNameB', ...
    'BESTALIGNED_STABILITY_THRESHOLD_CRITERIA', 'STABILITY_CLASSIFICATION_GROUP_LABEL', ...
    'BestAlignedStable', 'BestAlignedUnstable', ...
    'p1', 'p2', 'p3', 'p4', 'p5', 'p6');
fprintf('popvectorResults.mat saved to %s\n', outputFilename);

close(hFig)

end % function


%% FUNCTIONS
% close all
% 
% % This is for debugging and example figures
% T = get_maps(BestAlignedStable, 'Day 1', 'Day 3');
% maps1 = helper_map_cells_to_3dmatrix(T.mapDayAContext1);
% maps2 = helper_map_cells_to_3dmatrix(T.mapDayAContext2);
% [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
%             
% Plot some examples
% for iRow = 1:10 %size(T,1)
%     [hFig] = help_plot_T_row(T, iRow);
%     drawnow
% end % iRow


function [hFig] = help_plot_T_row(T, iRow)
hFig = figure()
subplot(2,2,1)
imagesc(T.mapDayAContext1{iRow})
title(sprintf('%s Context 1', T.sessionNameA{iRow}));
axis equal tight
set(gca, 'ydir', 'reverse')
colormap jet

subplot(2,2,2)
imagesc(T.mapDayAContext2{iRow})
title(sprintf('%s Context 2', T.sessionNameA{iRow}));
axis equal tight
set(gca, 'ydir', 'reverse')
colormap jet

subplot(2,2,3)
imagesc(T.mapDayBContext1{iRow})
title(sprintf('%s Context 1', T.sessionNameB{iRow}));
axis equal tight
set(gca, 'ydir', 'reverse')
colormap jet

subplot(2,2,4)
imagesc(T.mapDayBContext2{iRow})
title(sprintf('%s Context 2', T.sessionNameB{iRow}));
axis equal tight
set(gca, 'ydir', 'reverse')
colormap jet
sgtitle(sprintf('Registered Cell\n%s', T.registeredCellName{iRow}), 'interpreter', 'none');
end % function



% Compute the dot products
function [F] = compute_F(BestAligned)
    dayNames = {'Day 1', 'Day 2', 'Day 3'};
    numDays = length(dayNames);

    F = [];
    k = 1;
    % Compute for every unique pair of days.
    for iA = 1:numDays
        sessionNameA = dayNames{iA};
        for iB = iA+1:numDays
            sessionNameB = dayNames{iB};

            % Get the maps
            fprintf('Getting maps for %s and %s\n', sessionNameA, sessionNameB);
            T = get_maps(BestAligned, sessionNameA, sessionNameB);

            % Across A
            maps1 = helper_map_cells_to_3dmatrix(T.mapDayAContext1);
            maps2 = helper_map_cells_to_3dmatrix(T.mapDayAContext2);
            [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
            F(k).dp_across_day_A = dp;
            F(k).across_day_A_average = dpAvg;


            % Across B
            maps1 = helper_map_cells_to_3dmatrix(T.mapDayBContext1);
            maps2 = helper_map_cells_to_3dmatrix(T.mapDayBContext2);
            [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
            F(k).dp_across_day_B = dp;
            F(k).across_day_B_average = dpAvg;

            % Within C1
            maps1 = helper_map_cells_to_3dmatrix(T.mapDayAContext1);
            maps2 = helper_map_cells_to_3dmatrix(T.mapDayBContext1);
            [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
            F(k).dp_across_context_1 = dp;
            F(k).across_context_1_average = dpAvg;

            % Within C2
            maps1 = helper_map_cells_to_3dmatrix(T.mapDayAContext2);
            maps2 = helper_map_cells_to_3dmatrix(T.mapDayBContext2);
            [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
            F(k).dp_across_context_2 = dp;
            F(k).across_context_2_average = dpAvg;

    %         F(k).sessionNameA = sessionNameA;
    %         F(k).sessionNameB = sessionNameB;
            F(k).sessionNameA = sessionNameA;
            F(k).sessionNameB = sessionNameB;
            k = k + 1;
        end
    end
    F = struct2table(F);
end % function



function maps = helper_map_cells_to_3dmatrix(mapCells)
    m1 = mapCells{1};
    numMaps = length(mapCells);
    maps = nan(size(m1,1), size(m1,2), numMaps);
    for iMap = 1:numMaps
       m = mapCells{iMap};
       maps(:,:,iMap) = m(:,:); 
    end
end

function [T] = get_maps(BestAligned, sessionNameEarly, sessionNameLate)

    uniqueRegisteredCellNames = unique(BestAligned.registeredCellName); 
    numRegisteredCellNames = length(uniqueRegisteredCellNames);

    T = [];
    k = 1;
    for iName = 1:numRegisteredCellNames
        registeredCellName = uniqueRegisteredCellNames{iName};
        cdata = BestAligned(ismember(BestAligned.registeredCellName, registeredCellName) & ismember(BestAligned.groupLabel, {sessionNameEarly, sessionNameLate}),:);
        if size(cdata,1) ~= 2
            continue;
        end

        ind1 = find(ismember(cdata.groupLabel, sessionNameEarly));
        if length(ind1) ~= 1
            continue;
        end
        % Check that the map values are all valid
        if any(isnan(cdata.averageMapContext1{ind1}),'all') || any(isnan(cdata.averageMapContext2{ind1}),'all')
            continue;
        end
                
        ind2 = find(ismember(cdata.groupLabel, sessionNameLate));
        if length(ind2) ~= 1
            continue;
        end
        % Check that the map values are all valid
        if any(isnan(cdata.averageMapContext1{ind2}),'all') || any(isnan(cdata.averageMapContext2{ind2}),'all')
            continue;
        end
        
        
        T(k).registeredCellName = registeredCellName;
        T(k).mapDayAContext1 = cdata.averageMapContext1{ind1};
        T(k).mapDayAContext2 = cdata.averageMapContext2{ind1};
        T(k).mapDayBContext1 = cdata.averageMapContext1{ind2};
        T(k).mapDayBContext2 = cdata.averageMapContext2{ind2};
        T(k).sessionNameA = sessionNameEarly;
        T(k).sessionNameB = sessionNameLate;
        
        % We need to orient the maps to some standard
        % Rotate Day B context 1 to be most similar to Day A Context 1
%         mapA = T(k).mapDayAContext1;
%         mapB = T(k).mapDayBContext1;
%         mapBRot = rot90(mapB, 2);
%         c1 = corr(mapA(:), mapB(:));
%         c2 = corr(mapA(:), mapBRot(:));
%         if c2 > c1
%             % Then rotate BOTH day B maps
%             T(k).mapDayBContext1 = rot90(T(k).mapDayBContext1, 2);
%             T(k).mapDayBContext2 = rot90(T(k).mapDayBContext2, 2);
%         end

        k = k + 1;
    end
    T = struct2table(T);
end % function
    

%%
function [metric] = compute_max_0_180_correlation(mapA, mapB)
    mapBRot = rot90(mapB, 2);
    
    c1 = corr(mapA(:), mapB(:));
    c2 = corr(mapA(:), mapBRot(:));
    metric = max([c1, c2]);
end
