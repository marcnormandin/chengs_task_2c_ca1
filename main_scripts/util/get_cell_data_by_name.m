% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [output] = get_cell_data_by_name(T, animalName, sessionName, cellName)
    % Should be a table with one entry for each cells and each trial. eg.
    % if 2 cells and 12 trials each, table should have 24 rows.
    
    output = [];

    cellDataInds = find( ismember(T.animalName, animalName) & ismember(T.sessionName, sessionName) & ismember(T.cellName, cellName) );
    if isempty(cellDataInds)
        warning('Cell not found.')
    end

    possibleScalar3FieldNames = helper_possibleScalar3FieldNames();
    possibleScalar1FieldNames = helper_possibleScalar1FieldNames();
    for iField = 1:length(possibleScalar1FieldNames)
        fieldName = possibleScalar1FieldNames{iField};

        if ismember(fieldName, fieldnames(T))
            tmp = T.(fieldName)(cellDataInds);
            output.(sprintf('%ss',fieldName)) = reshape(tmp, length(tmp), 1);
        end
    end

    for iField = 1:length(possibleScalar3FieldNames)
        fieldName = possibleScalar3FieldNames{iField};

        if ismember(fieldName, fieldnames(T))
            xxx = T.(fieldName)(cellDataInds); % This will a cell array
            eg = xxx{1}; % example map to get the size
            yyy = nan(size(eg,1), size(eg,2), length(xxx));
            for i = 1:size(yyy,3)
                x = xxx{i};
                yyy(:,:,i) = x;
            end
            output.(sprintf('%ss',fieldName)) = yyy;
        end
    end

    if ismember('cellId', fieldnames(T))
        x = T.cellId(cellDataInds);
        if ~all(x)
            error('Cell ids are not all equal, but should be.');
        end
        output.cellId = x(1);
    end

    if ismember('cellName', fieldnames(T))
        cellNames = T.cellName(cellDataInds);
        for k = 2:length(cellNames)
            if ~strcmp(cellNames{k}, cellNames{1})
                error('Cell names are not equal, but should be.');
            end
        end
    end
    
    if ismember('dig', fieldnames(T))
       tmp = T.dig(cellDataInds);
       output.digs = reshape(tmp, length(tmp), 1);
    end
    
    animalSessionCellNames = T.animalSessionCellName(cellDataInds);
    if ~all(ismember(animalSessionCellNames, animalSessionCellNames{1}))
        error('Code bug present because all maps should be associated with the same animalSessionCellName, but they are not.\n');
    end
    animalSessionCellName = animalSessionCellNames{1}; % all will be the same

    output.animalName = animalName;
    output.sessionName = sessionName;
    output.cellName = cellName;
    output.animalSessionCellName = animalSessionCellName;
    output.cellDataInds = reshape(cellDataInds, length(cellDataInds), 1);
    output.numEntries = length(cellDataInds);
end % function
