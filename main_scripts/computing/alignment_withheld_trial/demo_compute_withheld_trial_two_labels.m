% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% Requires the MapsData table
MapsData = analysisInput.MapsData;

NORMALIZE_MAPS = true;
digTypesToCompare = {'Corr', 'Geo'};
clc

T = MapsData;
animalSessionCellNames = unique(T.animalSessionCellName);
nCells = length(animalSessionCellNames);

R = [];
for iCell = 1:nCells
    animalSessionCellName = animalSessionCellNames{iCell};

    Tc = T(ismember(T.animalSessionCellName, animalSessionCellName) & ismember(T.dig, digTypesToCompare),:);

    Maps = cat(3,Tc.map{:});
    mapLabels = Tc.dig;

    try
        withheldTable = compute_withheld_trial_two_labels(Maps, mapLabels, NORMALIZE_MAPS);
        
        withheldValid = withheldTable(withheldTable.isValid == true,:);
        numCorrect = sum(withheldValid.bestMatch == withheldValid.mapWithheldType);
        numPossible = size(withheldValid,1);
        ratioCorrect = numCorrect / numPossible;
        
        k = length(R) + 1;
        R(k).animalName = Tc.animalName{1};
        R(k).sessionName = Tc.sessionName{1};
        R(k).daynum = Tc.dayNum(1);
        R(k).cellName = Tc.cellName{1};
        R(k).animalSessionCellName = Tc.animalSessionCellName{1};
        R(k).numCorrect = numCorrect;
        R(k).numPossible = numPossible;
        R(k).ratioCorrect = ratioCorrect;
    catch e
        fprintf('Error processing %s\n', animalSessionCellName);
    end
end
R = struct2table(R)

%%
R = helper_add_group_labels_to_table(R, load_sessions_to_groups_table(analysisSettings), true)
groupLabels = {'Day 1', 'Day 2', 'Day 3'};
nGroups = length(groupLabels);
for iGroup = 1:nGroups
    gl = groupLabels{iGroup};
    Rg = R(ismember(R.groupLabel, gl),:);
    rMean(iGroup) = mean(Rg.ratioCorrect, 'omitnan');
end
figure
bar(1:nGroups, rMean);
xticks(1:nGroups);
xticklabels(groupLabels)


%%
close all
hFig = figure('position', get(0,'screensize'));


subplot(nMaps,5,k+1);
    plot_map(mapAverage0)
    title(sprintf('Average %s, (%0.2f)', uniqueLabels{1}, rho1))
    
    subplot(nMaps,5,k+2);
    plot_map(mapAverage1);
    title(sprintf('Average %s, (%0.2f)', uniqueLabels{2}, rho2))

    subplot(nMaps,5,k+3)
    plot_map(mapT)
    title(sprintf('Withheld %s', uniqueLabels{cT}))
    
    subplot(nMaps,5,k+4);
    scatter(mapAverage1(:), mapT(:), 5, 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'Markerfacealpha', 0.2);
    axis equal tight
    
    subplot(nMaps,5,k+5);
    scatter(mapAverage2(:), mapT(:), 5, 'b', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'Markerfacealpha', 0.2);
    axis equal tight
    
    k = k + 5;
    
    



    
function plot_map(M)
imagesc(M)
set(gca, 'ydir', 'reverse')
axis equal tight off
end
