% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R, success] = compute_rate_matrix_table(T, rateName, rateMatrixNumTrials)
    R = [];
    success = true;
    wb = [];
    
    % Get a table of unique sets of (animal name, session name, cell name)
    records = get_sessions_and_cellnames_from_table(T);
    numRecords = size(records,1);

    wb = waitbar(0, 'Starting', 'name', sprintf('Computing Rate Matrices (%s)', rateName));
        
    k = 1;
    
    fprintf('Processing (%d) Rate Matrices (%s): ', numRecords, rateName);
    for iRec = 1:numRecords
        waitbar(iRec/numRecords, wb, sprintf('Progress: %d %%', floor(iRec/numRecords*100)))

        
        animalName = records.animalName{iRec};
        sessionName = records.sessionName{iRec};
        cellName = records.cellName{iRec};
       
        
        cellData = get_cell_data_by_name(T, animalName, sessionName, cellName);

        cellRateMatrixData = compute_cell_rate_difference_matrix_by_context(cellData.(rateName), cellData.trialIds, cellData.contextIds, rateMatrixNumTrials);

        try
            R(k).animalName = animalName;
            R(k).sessionName = sessionName;

            R(k).cellName = cellName;
            R(k).cellId = cellData.cellId;
            
            R(k).animalSessionCellName = cellData.animalSessionCellName;

            R(k).numRates = cellRateMatrixData.numRates;
            R(k).rateName = rateName;
            R(k).rateMatrix = cellRateMatrixData.rateMatrix;
            R(k).rateMatrixTrialIds = cellRateMatrixData.rateMatrixTrialIds;
            R(k).rateMatrixContextIds = cellRateMatrixData.rateMatrixContextIds;

            R(k).numRatesContext1 = sum(cellData.contextIds == 1);
            R(k).numRatesContext2 = sum(cellData.contextIds == 2);

            k = k + 1;
        catch e
            success = false;
            fprintf(1, 'Error while processing %s %s %s: %s\n', animalName, sessionName, cellName, e.identifier);
            fprintf(1, '\t%s\n', e.message);
        end
    end
    R = struct2table(R);
    
    fprintf(' done!\n');
    
    if ~isempty(wb)
        close(wb);
    end
end % function
