% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [output] = ml_alg_best_aligned_average_context_maps_for_cell_v2(maps, trialIds, contextIds, comparisonMethod, REFLECT_ONE_CONTEXT)
%% 
% sdataMaps should be a MxNxK matrix where K is the number of cells
% sdataCellIds should be a Kx1 array where each element is the cell id
% corresponding to the map
% sdataTrialIds should be a Kx1 array where each element is the trial id
% for the corresponding map
% sdataContextIds should be a Kx1 array where each elemen is the context id
% for the corresponding map

    numMaps = size(maps,3);

    for iMap = 1:numMaps
       mm = maps(:,:,iMap);
       mm(~isfinite(mm)) = 0;

       maps(:,:,iMap) = mm;
    end

    cellMaps1 = maps(:,:, contextIds == 1);
    cellMaps2 = maps(:,:, contextIds == 2);
    
    [finalMean1, bestCombination1, metricValue1] = best_flip_combinations(cellMaps1, comparisonMethod);
    [finalMean2, bestCombination2, metricValue2] = best_flip_combinations(cellMaps2, comparisonMethod);
    
    if REFLECT_ONE_CONTEXT
        bestCombination2 = flip(bestCombination2,2);
    end

    c1 = corr(finalMean1(:), finalMean2(:));
    c2 = corr(finalMean1(:), rot90(finalMean2(:), 2));
    finalCorrelation = c1;
    if c2 > c1
        finalMean2 = rot90( finalMean2, 2 );
        bestCombination2 = xor(bestCombination2, ones(1,length(bestCombination2)));
        finalCorrelation = c2;
    end

    flippedMaps1 = cellMaps1;
    ind1 = find(bestCombination1 == 1);
    for k = 1:length(ind1)
        flippedMaps1(:,:,ind1(k)) = rot90( flippedMaps1(:,:,ind1(k)), 2 );
    end
    flippedMaps2 = cellMaps2;
    ind2 = find(bestCombination2 == 1);
    for k = 1:length(ind2)
        flippedMaps2(:,:,ind2(k)) = rot90( flippedMaps2(:,:,ind2(k)), 2 );
    end
    
    % store the output
    output.context1.meanMap = finalMean1;
    output.context1.flipSequence = double(bestCombination1);
    output.context1.flippedMaps = flippedMaps1;
    output.context1.trialIds = trialIds( contextIds == 1 );
    output.context1.metricValue = metricValue1;
    
    output.context2.meanMap = finalMean2;
    output.context2.flipSequence = double(bestCombination2);
    output.context2.flippedMaps = flippedMaps2;
    output.context2.trialIds = trialIds( contextIds == 2 );
    output.context2.metricValue = metricValue2;

    output.bestCorrelation = finalCorrelation;

end % function


