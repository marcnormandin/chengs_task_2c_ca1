% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024

function [predictionDataset] = filter_prediction_dataset(predictionDataset, params)
    nRows = size(predictionDataset,1);
    nCellsBeforeFiltering = nan(nRows,1);
    nCellsAfterFiltering = nan(nRows,1);
    emptyRows = [];
    for iRow = 1:nRows
        
        bestCorrelations = predictionDataset.bestCorrelations{iRow};

        indsValid = 1:length(bestCorrelations); % all are initially valid

        nCellsBeforeFiltering(iRow) = length(bestCorrelations);
        if strcmpi(params.stabilityType, 'stable')
            indsValidStabilityType = find(bestCorrelations >= params.stabilityThreshhold);
        elseif strcmpi(params.stabilityType, 'unstable')
            indsValidStabilityType = find(bestCorrelations < params.stabilityThreshhold);
        elseif strcmpi(params.stabilityType, 'any')
            indsValidStabilityType = 1:length(bestCorrelations); % use all
        else
            error('stabilityType must be stable, unstable or any');
        end
        indsValid = intersect(indsValid, indsValidStabilityType);

        IS_CALCIUM_DATA = any(ismember(predictionDataset.Properties.VariableNames, {'cellRegScores'}));
        if IS_CALCIUM_DATA
            cellRegScores = predictionDataset.cellRegScores{iRow};
            indsValidCellRegScore = find(cellRegScores >= params.minCellRegScore);
            indsValid = intersect(indsValid, indsValidCellRegScore);
        end

        % Finding valid indices is now done. Now use them to filter.

        
        if IS_CALCIUM_DATA
            cellRegScores = predictionDataset.cellRegScores{iRow};
            cellRegScores = cellRegScores(indsValid);
            predictionDataset.cellRegScores{iRow} = cellRegScores; % store
        end

        maps = predictionDataset.maps{iRow};
        cellNames = predictionDataset.cellNames{iRow};
        centerOutAnglesDeg = predictionDataset.centerOutAnglesDeg{iRow};
        % bestCorrelations already retrieved
        regCellNames = predictionDataset.regCellNames{iRow};
        radialDistance = predictionDataset.radialDistance{iRow};

        % Reduce
        maps = maps(:,:,indsValid);
        cellNames = cellNames(indsValid);
        centerOutAnglesDeg = centerOutAnglesDeg(indsValid);
        bestCorrelations = bestCorrelations(indsValid);
        regCellNames = regCellNames(indsValid);
        radialDistance = radialDistance(indsValid);

        nCellsAfterFiltering(iRow) = length(bestCorrelations);

        % Store
        predictionDataset.maps{iRow} = maps;
        predictionDataset.cellNames{iRow} = cellNames;
        predictionDataset.centerOutAnglesDeg{iRow} = centerOutAnglesDeg;
        predictionDataset.bestCorrelations{iRow} = bestCorrelations;
        predictionDataset.regCellNames{iRow} = regCellNames;
        predictionDataset.radialDistance{iRow} = radialDistance;

        if nCellsAfterFiltering(iRow) == 0 || nCellsAfterFiltering(iRow) <= params.minCells
            emptyRows(end+1) = iRow;
        end
    end % iRow
    predictionDataset.nCellsBeforeFiltering = nCellsBeforeFiltering;
    predictionDataset.nCellsAfterFiltering = nCellsAfterFiltering;

    % If no cells remain, then remove the row.
    predictionDataset(emptyRows,:) = [];
end % function
