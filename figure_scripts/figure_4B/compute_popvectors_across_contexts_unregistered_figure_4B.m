% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function compute_popvectors_across_contexts_unregistered_figure_4B(analysisSettings, BestAligned, BESTALIGNED_STABILITY_THRESHOLD_CRITERIA, OUTPUT_FOLDER)
    % Add the group labels (day names) to the data
    SessionsToGroups = load_sessions_to_groups_table(analysisSettings);
    eliminateUnmatched = true;
    BestAligned = helper_add_group_labels_to_table(BestAligned, SessionsToGroups, eliminateUnmatched);
    fprintf('Added the day/group label information to the BestAligned results.\n');
    
    groupLabels = {'Day 1', 'Day 2', 'Day 3'};
    
    % We only care about these.
    BestAligned(~ismember(BestAligned.groupLabel, groupLabels),:) = [];
    
    BestAlignedStable = BestAligned(BestAligned.bestCorrelation >= BESTALIGNED_STABILITY_THRESHOLD_CRITERIA,:);
    BestAlignedUnstable = BestAligned(BestAligned.bestCorrelation < BESTALIGNED_STABILITY_THRESHOLD_CRITERIA,:);
    
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
    sessionNamesToShow = {'Day 1', 'Day 2', 'Day 3'};
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
    xlabel('Dot Product')
    ylabel('Cumulative')
    
    FSz = 24;
    
    xticks(0:0.2:1)
    yticks(0:0.2:1);
    set(gca, 'fontsize', FSz)
    axis square
    
    % save the figure
    mulana_savefig(hFig, OUTPUT_FOLDER, sprintf('figure_4B_calcium_popvectors_unregistered_across_context_classified_per_day'), {'png', 'svg', 'fig'});
    
    
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
    
    writetable(FS, fullfile(OUTPUT_FOLDER, 'figure_4B_calcium_popvectors_unregistered_across_context_classified_per_day.xlsx'));
    save(fullfile(OUTPUT_FOLDER, 'figure_4B_calcium_popvectors_unregistered_across_context_classified_per_day.mat'), 'analysisSettings', 'F', 'FS')
end % main function



%% HELPER FUNCTIONS

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
