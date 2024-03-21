% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [F, success] = ml_alg_compute_best_aligned_map_table_v3(T, comparisonMethod, REFLECT_ONE_CONTEXT)
    F = [];
    success = true;
    wb = [];
    
    % Get unique animalName, sessionName, cellName combinations
    R = get_sessions_and_cellnames_from_table(T);
    numRecords = size(R,1);

    fprintf('Processing (%d) Best Aligned Maps and correlations: ', numRecords);

    wb = waitbar(0, 'Starting', 'name', 'Computing Best Aligned Maps and correlations');
        
    k = 1;
    % Unique cell record
    for iRec = 1:numRecords
        waitbar(iRec/numRecords, wb, sprintf('Progress: %d %%', floor(iRec/numRecords*100)))
        
        animalName = R.animalName{iRec};
        sessionName= R.sessionName{iRec};
        dayNum = R.dayNum(iRec);
        cellName = R.cellName{iRec};
        
        try
            cellData = get_cell_data_by_name(T, animalName, sessionName, cellName);
            
            % If we are reflecting a contexts maps
%             if ~isempty(BESTALIGNED_CONTEXTS_TO_REFLECT_MAPS)
%                 for iContext = 1:length(BESTALIGNED_CONTEXTS_TO_REFLECT_MAPS)
%                     contextId = BESTALIGNED_CONTEXTS_TO_REFLECT_MAPS(iContext);
%                     matchingInds = find(cellData.contextIds == contextId);
%                     for j = 1:length(matchingInds)
%                         m = cellData.maps(:,:,matchingInds(j));
%                         cellData.maps(:,:,matchingInds(j)) = flip(m,2);
%                     end
%                 end
%             end
            
            
            output = ml_alg_best_aligned_average_context_maps_for_cell_v3(cellData.maps, cellData.trialIds, cellData.contextIds, comparisonMethod, REFLECT_ONE_CONTEXT);

            F(k).animalName = animalName;
            F(k).sessionName = sessionName;
            F(k).dayNum = dayNum;
            F(k).cellName = cellName;
            F(k).animalSessionCellName = cellData.animalSessionCellName;
            F(k).cellId = cellData.cellId;
            F(k).comparisonMethod = comparisonMethod;
            F(k).cellDataInds = cellData.cellDataInds;
            
            F(k).metricValue1 = output.context1.metricValue;
            F(k).metricValue2 = output.context2.metricValue;
            F(k).metricValueAll = output.all.metricValue;
            
            F(k).averageMapContext1(:,:) = output.context1.meanMap;
            F(k).averageMapContext2(:,:) = output.context2.meanMap;
            F(k).averageMapContextAll(:,:) = output.all.meanMap;
            
            F(k).rotationSequence1 = output.context1.flipSequence;
            F(k).rotationSequence2 = output.context2.flipSequence;
            F(k).rotationSequenceAll = output.all.flipSequence;
            
            F(k).trialIdsContext1 = output.context1.trialIds;
            F(k).trialIdsContext2 = output.context2.trialIds;

            F(k).numTrialsContext1 = length(output.context1.trialIds);
            F(k).numTrialsContext2 = length(output.context2.trialIds);
            F(k).numTrialsContextAll = length(output.all.trialIds);

            F(k).maxTrials = max([length(output.context1.trialIds), length(output.context2.trialIds)]);

            F(k).bestCorrelation = output.bestCorrelation;
                        
            k = k + 1;
        catch e
            success = false;
            fprintf(1, 'Error while processing %s %s %s: %s\n', animalName, sessionName, cellName, e.identifier);
            fprintf(1, '\t%s\n', e.message);
        end
    end

    F = struct2table(F);
    
    fprintf(' done!\n');
    
    if ~isempty(wb)
        close(wb);
    end
end % function

