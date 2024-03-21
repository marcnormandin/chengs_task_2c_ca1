% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [averageCorrelations, totalSamples, numCells, numSamplesExcluded, animalSessionCellIds, animalSessionCellNames] = compute_bfo180_similarity_using_percell_table(PerCell180, fieldNamesWanted, BFO180_MIN_MAPS)
    BFO180_MIN_CORRELATION = []; % 
    
    fieldNamesPossible = {'any', 'context1', 'context2', 'different'};
    if ~all(ismember(fieldNamesWanted, fieldNamesPossible))
        error('Invalid field name(s) requested.\n');
    end
    numFieldNames = length(fieldNamesWanted);
    
    numCells = size(PerCell180,1);

    % Compute minimum number of comparisons based on the minimum number of
    % trials.
    MIN_COMPARISONS = BFO180_MIN_MAPS*(BFO180_MIN_MAPS-1)/2;
    
    % Store the best found angle indices
    best_inds = [];
    best_v = [];
    
    not_used_vinds = []; % equivalent of best_inds except we don't use them because of the inclusion criteria
    not_used_vs = [];
    
    % Each cell will contribute a single average correlation value
    averageCorrelations = [];
    animalSessionCellIds = [];
    animalSessionCellNames = {};
    for iCell = 1:numCells
       vind = [];
       v = [];
       asid = PerCell180.animalSessionCellId(iCell);
       ascn = PerCell180.animalSessionCellName{iCell};
       
       % This allow for combining c1-c1 and c2-c2 results to get 'within'
       for iFieldName = 1:numFieldNames
            fieldName = fieldNamesWanted{iFieldName};
            
            % The angle index
            vind_cell = PerCell180.(sprintf('vind_%s', fieldName)){iCell};
            % The correlation for the angle
            v_cell = PerCell180.(sprintf('v_%s', fieldName)){iCell};
            
            vind_cell = reshape(vind_cell, length(vind_cell), 1);
            v_cell = reshape(v_cell, length(v_cell), 1);
            
            badInd = find(~isfinite(vind_cell) | ~isfinite(v_cell));
            not_used_vinds = cat(1, not_used_vinds, vind_cell(badInd));
            not_used_vs = cat(1, not_used_vs, v_cell(badInd));
            vind_cell(badInd) = [];
            v_cell(badInd) = [];
            
            if ~isempty(BFO180_MIN_CORRELATION)
                vind_use = vind_cell(v_cell >= BFO180_MIN_CORRELATION);
                v_use = v_cell(v_cell >= BFO180_MIN_CORRELATION);
                
                not_used_vinds = cat(1, not_used_vinds, vind_cell(v_cell < BFO180_MIN_CORRELATION));
                not_used_vs = cat(1, not_used_vs, v_cell(v_cell < BFO180_MIN_CORRELATION));
            else
                vind_use = vind_cell; % use all the data
                v_use = v_cell;
            end
            
            % The correlation value for the best angle found
            %vind = [vind, vind_use];
            %v = [v, v_use];
            vind = cat(1, vind, vind_use);
            v = cat(1, v, v_use);
       end
       
       % Only use data if it meets the minimum number of comparisons
       % criteria
       if length(vind) < MIN_COMPARISONS
           not_used_vs = cat(1, not_used_vs, v);
           not_used_vinds = cat(1, not_used_vinds, vind);
           continue;
       end
       
       %best_inds = [best_inds, vind];
       %best_v = [best_v, v];
       best_inds = cat(1, best_inds, vind);
       best_v = cat(1, best_v, v);
       averageCorrelations = cat(1, averageCorrelations, nanmean(v));
       animalSessionCellIds = cat(1, animalSessionCellIds, asid);
       animalSessionCellNames = cat(1, animalSessionCellNames, ascn);
    end
    
    %averageCorrelations = best_v;
    
    totalSamples = length(best_v);
    
    % not used
    numSamplesExcluded = length(not_used_vinds);
end
