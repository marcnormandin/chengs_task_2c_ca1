% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [output] = ml_algo_bfo_percell_general_single(maps, contextIds, rotationsDeg, MIN_DIFFERENCE_FOR_MAX_CORRELATION)
    % Force to use two contexts because some cells (almost none) only have data for one.
    uniqueContextIds = sort(unique(contextIds));
    numContexts = length(uniqueContextIds);
    
    % Check that we can apply this function to the maps
    examplePlacemap = maps(:,:,1);
    placemapDim1 = size(examplePlacemap,1);
    placemapDim2 = size(examplePlacemap,2);
    placemapIsSquare = (placemapDim1 == placemapDim2);
    if any(ismember(rotationsDeg, [90, 270])) && ~placemapIsSquare
        error('Can not proceed because 90 or 270 degree rotations are requested, BUT the maps are not square. Maps are (%d, %d)', placemapDim1, placemapDim2);
    end

    % Now separate maps based on contexts.
    contextMaps = cell(numContexts,1);
    for iContext = 1:numContexts
        contextMaps{iContext} = maps(:,:, contextIds == uniqueContextIds(iContext));
    end

    % Any two contexts (all), c1-c1, c1-c2, c2-c2
    [output.v_any, output.vind_any, rotationsDeg] = ml_alg_bfo(rotationsDeg, maps, maps, false, MIN_DIFFERENCE_FOR_MAX_CORRELATION);

    % Same context, eg c1-c1, c2-c2
    for iContext = 1:numContexts
        [output.(sprintf('v_context%d', uniqueContextIds(iContext))), output.(sprintf('vind_context%d', uniqueContextIds(iContext))), rotationsDeg] ...
            = ml_alg_bfo(rotationsDeg, contextMaps{iContext}, contextMaps{iContext}, false, MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    end

    % 
    if numContexts ~= 2
        warning('This only works for two contexts')
        output.v_different = [];
        output.vind_different = [];
    else
        [output.v_different, output.vind_different, rotationsDeg] = ml_alg_bfo(rotationsDeg, contextMaps{1}, contextMaps{2}, true, MIN_DIFFERENCE_FOR_MAX_CORRELATION);
    end
    
end % function
