% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function compute_popvectors_across_contexts_unregistered(analysisSettings, BestAligned, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA)


% Add the group labels (day names) to the data
SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
eliminateUnmatched = true;
BestAligned = helper_add_group_labels_to_table(BestAligned, SessionsToGroups, eliminateUnmatched);
fprintf('Added the day/group label information to the BestAligned results.\n');

% TEMP
% BestAligned(BestAligned.metricValue1 <= 0.5 | BestAligned.metricValue2 <= 0.5,:) = [];
% warning('using a temp, not final!!!');

groupLabels = {'Day 1', 'Day 2', 'Day 3'};

% % Classify the stability of the cells based on the given day.
% [registeredCellNamesStable, registeredCellNamesUnstable] =  get_stable_unstable_registered_cell_names_classified_1day(analysisSettings, CellRegTable, BestAligned, STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);
% fprintf('There are %d stable cells, and %d unstable cells, based on %s with a minimum correlation threshold of %0.1f\n', length(registeredCellNamesStable), length(registeredCellNamesUnstable), STABILITY_CLASSIFICATION_GROUP_LABEL, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA);
% 
% % Add the associated registered cell names to the best aligned table
% BestAligned = bestaligned_add_registered_cell_names(BestAligned, CellRegTable);

% We only care about these "Days".
BestAligned(~ismember(BestAligned.groupLabel, groupLabels),:) = [];

% Note, that we might want to update the 'isStable' column of the
% BestAligned table, but the code works without it because we select by
% stable and unstable registered cell names and that used the input
% correlation from the main function input, not that necessarily set by the
% main analysis settings.
%BestAlignedStable = BestAligned(ismember(BestAligned.registeredCellName, registeredCellNamesStable),:);
%BestAlignedUnstable = BestAligned(ismember(BestAligned.registeredCellName, registeredCellNamesUnstable), :);

BestAlignedStable = BestAligned(BestAligned.bestCorrelation >= BESTALIGNED_STABILITY_THRESHOLD_CRITERIA,:);
BestAlignedUnstable = BestAligned(BestAligned.bestCorrelation < BESTALIGNED_STABILITY_THRESHOLD_CRITERIA,:);


% % Go through each day and compute the cell type counts
% for iGroup = 1:length(groupLabels)
%     groupLabel = groupLabels{iGroup};
%     T1 = BestAlignedStable(ismember(BestAlignedStable.groupLabel, groupLabel),:);
%     T2 = BestAlignedUnstable(ismember(BestAlignedUnstable.groupLabel, groupLabel),:);
%     numStable = size(T1,1);
%     numUnstable = size(T2,1);
%     numTotal = numStable + numUnstable;
%     fprintf('Stable cells make up %0.2f of cells on %s\n', numStable/numTotal*100, groupLabel);
% end



% Compute stuff with bestaligned maps across days

% Compute the within and across for all combinations of days for both types
% of cells.
[FStable] = compute_F(BestAlignedStable);
[FUnstable] = compute_F(BestAlignedUnstable);

FStable.isStable = true(size(FStable,1),1);
FUnstable.isStable = false(size(FUnstable,1),1);
%%
F = [FStable; FUnstable];
F.cumdist_x = cell(size(F,1),1);
F.cumdist_y = cell(size(F,1),1);
for iRow = 1:size(F,1)
    [x,y] = ml_alg_cumdist(F.dp_across{iRow});
    F.cumdist_x{iRow} = x;
    F.cumdist_y{iRow} = y;
end % iRow

close all

hFig = figure('position', get(0,'screensize'));
sessionNames = {'Day 1', 'Day 2', 'Day 3'};
sessionNamesToShow = {'Day 1', 'Day 3'};
sessionColours = {'r','g','b'};
stabilityLineStyles = {':','-'};

LW = 4;
l = {};
for iRow = 1:size(F,1)
    sessionName = F.sessionName{iRow};
    isStable = F.isStable(iRow);
    
    if ~ismember(sessionNamesToShow, sessionName)
        continue;
    end
    plot(F.cumdist_x{iRow}, F.cumdist_y{iRow}, 'color', sessionColours{find(ismember(sessionNames, sessionName))}, 'LineStyle', stabilityLineStyles{isStable+1}, 'linewidth', LW);
    hold on
    if isStable
        isStableStr = 'stable';
    else
        isStableStr = 'unstable';
    end
    l{end+1} = sprintf('%s across %s', sessionName, isStableStr);
end
legend(l, 'location', 'southoutside')
grid on
grid minor
%sgtitle(sprintf('Population Vectors Analysis\nAll cells (non-registered'))
xlabel('Dot Product')
ylabel('Cumulative')

FSz = 24;

xticks(0:0.2:1)
yticks(0:0.2:1);
set(gca, 'fontsize', FSz)
axis square


% save the output to here
outputFolder = fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'popvectors');
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% save the figure
mulana_savefig(hFig, outputFolder, sprintf('calcium_popvectors_unregistered_across_context_classified_per_day'), {'png', 'svg', 'fig'});

% Save the data
% outputFilename = fullfile(outputFolder, sprintf('calciumPopvectorResults_unregistered_across_context_classified_per_day.mat') );
% save(outputFilename, 'popvectorResults', 'analysisSettings', ...
%     'sessionNameA', 'sessionNameB', ...
%     'BESTALIGNED_STABILITY_THRESHOLD_CRITERIA', 'STABILITY_CLASSIFICATION_GROUP_LABEL', ...
%     'BestAlignedStable', 'BestAlignedUnstable', ...
%     'p1', 'p2', 'p3', 'p4', 'p5', 'p6');
% fprintf('popvectorResults.mat saved to %s\n', outputFilename);

%close(hFig)





% These compare stable with stable and unstable with unstable
% [ks1,p1] = kstest2(FS.('dp_across_day_A'){1},FS.('dp_across_day_B'){1});
% [ks2,p2] = kstest2(FU.('dp_across_day_A'){1},FU.('dp_across_day_B'){1});
% 
% % These compare stable with unstable
% [ks3,p3] = kstest2(FS.('dp_across_day_A'){1},FU.('dp_across_day_A'){1});
% [ks4,p4] = kstest2(FS.('dp_across_day_B'){1},FU.('dp_across_day_B'){1});
% [ks5,p5] = kstest2(FS.('dp_across_day_A'){1},FU.('dp_across_day_B'){1});
% [ks6,p6] = kstest2(FS.('dp_across_day_B'){1},FU.('dp_across_day_A'){1});


%% KS-tests
FS = [];
k = 1;
for iRow1 = 1:size(F,1)
    x1 = F.('dp_across'){iRow1};
    for iRow2 = iRow1+1:size(F,1)
        x2 = F.('dp_across'){iRow2};
        [ks1,p1] = kstest2(x1,x2); 
        FS(k).kstest_H0 = ks1;
        FS(k).kstest_p = p1;
        
        FS(k).numAnimalsA = F.numAnimals(iRow1);
        FS(k).numCellsA = F.numCells(iRow1);
        FS(k).sessionNameA = F.sessionName{iRow1};
        FS(k).isStableA = F.isStable(iRow1);
        
        FS(k).numAnimalsB = F.numAnimals(iRow2);
        FS(k).numCellsB = F.numCells(iRow2);
        FS(k).sessionNameB = F.sessionName{iRow2};
        FS(k).isStableB = F.isStable(iRow2);
        
        k = k + 1;
    end
    
end
FS = struct2table(FS);

writetable(FS, fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'popvectors', 'calcium_popvectors_unregistered_across_context_classified_per_day.xlsx'));

save(fullfile(analysisSettings.OUTPUT_RESULTS_FOLDER, 'popvectors', 'calcium_popvectors_unregistered_across_context_classified_per_day.mat'), 'analysisSettings', 'F', 'FS')

end % main function



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
        sessionName = dayNames{iA};
        
        % Get the maps
        fprintf('Getting maps for %s\n', sessionName);
        T = get_maps(BestAligned, sessionName);

        % Across A
        maps1 = helper_map_cells_to_3dmatrix(T.mapContext1);
        maps2 = helper_map_cells_to_3dmatrix(T.mapContext2);
        [dp, dpAvg] = ml_alg_popvectors_compute_dotproducts(maps1, maps2);
        F(k).dp_across = dp;
        F(k).dp_across_average = dpAvg;
        F(k).numCells = size(maps1,3); % same as maps2
        F(k).numAnimals = length(unique(T.animalName));

        F(k).sessionName = sessionName;
        k = k + 1;
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

% Get all the maps that match the session, that have a context 1 and
% context 2 map.
function [T] = get_maps(BestAligned, sessionName)

    uniqueCellNames = unique(BestAligned.animalSessionCellName); 
    numCellNames = length(uniqueCellNames);

    T = [];
    k = 1;
    for iName = 1:numCellNames
        cellName = uniqueCellNames{iName};
        cdata = BestAligned(ismember(BestAligned.animalSessionCellName, cellName) & ismember(BestAligned.groupLabel, sessionName),:);
        if size(cdata,1) ~= 1
            continue;
        end

        ind1 = find(ismember(cdata.groupLabel, sessionName));
        if length(ind1) ~= 1
            continue;
        end
        % Check that the map values are all valid
        if any(isnan(cdata.averageMapContext1{ind1}),'all') || any(isnan(cdata.averageMapContext2{ind1}),'all')
            continue;
        end
                
        ind2 = find(ismember(cdata.groupLabel, sessionName));
        if length(ind2) ~= 1
            continue;
        end
        % Check that the map values are all valid
        if any(isnan(cdata.averageMapContext1{ind2}),'all') || any(isnan(cdata.averageMapContext2{ind2}),'all')
            continue;
        end
        
        
        T(k).animalSessionCellName = cellName;
        T(k).mapContext1 = cdata.averageMapContext1{ind1};
        T(k).mapContext2 = cdata.averageMapContext2{ind1};
        T(k).sessionName = sessionName;
        T(k).animalName = cdata.animalName{1};
        
        % We need to orient the maps to some standard
        % Rotate Day B context 1 to be most similar to Day A Context 1
        map1 = T(k).mapContext1;
        map2 = T(k).mapContext2;
        map2Rot = rot90(map2, 2);
        c1 = corr(map1(:), map2(:));
        c2 = corr(map1(:), map2Rot(:));
        if c2 > c1
            % Then rotate BOTH day B maps
            T(k).mapContext2 = map2Rot;
        end

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
