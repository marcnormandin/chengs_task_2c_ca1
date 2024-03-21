% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [PerCell, success] = compute_percell_for_table(MapsData, ROTATIONS_DEG, CONTEXTS_TO_REFLECT_MAPS, MIN_DIFFERENCE_FOR_MAX_CORRELATION)
    PerCell = [];
    success = true;
    wb = [];
    
    % ROTATIONS_DEG = [0, 90, 180, 270] or [0, 180]
    SessionsAndCellsTable = get_sessions_and_cellnames_from_table(MapsData);
    
    numRecords = size(SessionsAndCellsTable,1);
    
    fprintf('Processing (%d) Per Cell records: ', numRecords);
    wb = waitbar(0, 'Starting', 'name', 'Computing PerCell');
    
    k = 1;
    for iRec = 1:numRecords
        waitbar(iRec/numRecords, wb, sprintf('Progress: %d %%', floor(iRec/numRecords*100)))
        
        animalName = SessionsAndCellsTable.animalName{iRec};
        sessionName = SessionsAndCellsTable.sessionName{iRec};
        cellName = SessionsAndCellsTable.cellName{iRec};
        
%         try
            cellData = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
            
            % If we are reflecting a contexts maps
            if ~isempty(CONTEXTS_TO_REFLECT_MAPS)
                for iContext = 1:length(CONTEXTS_TO_REFLECT_MAPS)
                    contextId = CONTEXTS_TO_REFLECT_MAPS(iContext);
                    matchingInds = find(cellData.contextIds == contextId);
                    fprintf('\t\treflecting context %d\n', contextId);
                    for j = 1:length(matchingInds)
                        m = cellData.maps(:,:,matchingInds(j));
                        cellData.maps(:,:,matchingInds(j)) = fliplr(m); % this was flip(m,2) which is the same
                    end
                end
            end

            perCell = ml_algo_bfo_percell_general_single(cellData.maps, cellData.contextIds, ROTATIONS_DEG, MIN_DIFFERENCE_FOR_MAX_CORRELATION);

            PerCell(k).animalName = animalName;
            PerCell(k).sessionName = sessionName;
            PerCell(k).cellName = cellName;
            PerCell(k).cellId = cellData.cellId;
            PerCell(k).rotations = ROTATIONS_DEG;
            PerCell(k).animalSessionCellName = cellData.animalSessionCellName;
            
            fieldNames = fieldnames(perCell);
            for iField = 1:length(fieldNames)
                fieldName = fieldNames{iField};
                PerCell(k).(fieldName) = perCell.(fieldName);
            end
            
            k = k + 1;
%         catch e
%             success = false;
%             fprintf(1, 'Error while processing %s %s %s: %s\n', animalName, sessionName, cellName, e.identifier);
%             fprintf(1, '\t%s\n', e.message);
%         end
    end
    PerCell = struct2table(PerCell);
    
    fprintf(' done!\n');
    
    if ~isempty(wb)
        close(wb);
    end
end % function
