% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [analysisInput] = analysis_load_marc_calcium_input_data(analysisSettings, dataFullFilename)
    % Load the main rectangular placemap data
    if isempty(dataFullFilename)
        [dataFilename, dataPath] = uigetfile();
        dataFullFilename = fullfile(dataPath, dataFilename);
    end
    tmp = load(dataFullFilename);
    
    MapsData = tmp.MapsData;
    MapsSquareData = tmp.MapsSquareData;
    
    % Add the unique cell names
    uniqueCellNames = cell(size(MapsData,1),1);
    for i = 1:size(MapsData,1)
        animalName = MapsData.animalName{i};
        sessionName = MapsData.sessionName{i};
        cellName = MapsData.cellName{i};
        uniqueCellNames{i} = sprintf('%s_%s_%s', animalName, sessionName, cellName);
    end
    MapsData.animalSessionCellName = uniqueCellNames;
    
    uniqueCellNames = cell(size(MapsSquareData,1),1);
    for i = 1:size(MapsSquareData,1)
        animalName = MapsSquareData.animalName{i}
        sessionName = MapsSquareData.sessionName{i}
        cellName = MapsSquareData.cellName{i}
        uniqueCellNames{i} = sprintf('%s_%s_%s', animalName, sessionName, cellName);
    end
    MapsSquareData.animalSessionCellName = uniqueCellNames;
    
    
    % If needed, use the settings to exclude maps
    MapsData.map = MapsData.(analysisSettings.CALCIUM_MAPNAME_TO_USE);
    MapsSquareData.map = MapsSquareData.(analysisSettings.CALCIUM_MAPNAME_TO_USE);
    
    % If set, eliminate data where animal did not dig
    if analysisSettings.ELIMINATE_NO_DIGS
       inds = find(~ismember(MapsData.dig, {'Corr', 'Geo', 'Feat', 'Wrong'}));
       MapsData(inds,:) = [];

       inds = find(~ismember(MapsSquareData.dig, {'Corr', 'Geo', 'Feat', 'Wrong'}));
       MapsSquareData(inds,:) = [];
    end
    clear inds
    
    % Use the filters to eliminate bad data
    if analysisSettings.CALCIUM_APPLY_CELLLIKE_SPATIAL_FOOTPRINT_FILTER
        inds = find(MapsData.passedCelllikeSpatialFootprintFilter == 0);
        MapsData(inds,:) = [];
        
        inds = find(MapsSquareData.passedCelllikeSpatialFootprintFilter == 0);
        MapsSquareData(inds,:) = [];
    end
    
    MapsDataExcluded = [];
    MapsSquareDataExcluded = [];
    
    % Filter, if non-empty, using the information content
    if ~isempty(analysisSettings.CALCIUM_ICS_FILTER_PERCENTILE)
        fprintf('Filtering out all maps whose information content (ICS) is lower than the %d-th percentile.\n', analysisSettings.CALCIUM_ICS_FILTER_PERCENTILE);
        % Rectangular maps
        ics = real(MapsData.ics(:));
        icsV = prctile(ics, analysisSettings.CALCIUM_ICS_FILTER_PERCENTILE);
        
        MapsDataExcluded = MapsData(ics < icsV,:);
        
        MapsData(ics < icsV,:) = [];
        
       
        % Rectangular maps
        ics = real(MapsSquareData.ics(:));
        icsV = prctile(ics, analysisSettings.CALCIUM_ICS_FILTER_PERCENTILE);
        
        MapsSquareDataExcluded = MapsSquareData(ics < icsV,:);
        MapsSquareData(ics < icsV,:) = [];
    else
        fprintf('No filtering based on information content applied.\n');
    end
        
    
    % Add the cellreg data
    [dataCellRegFilename, dataCellRegPath] = uigetfile({'*.mat'}, 'Select calcium_CellRegTable mat file.');
    dataCellRegFullFilename = fullfile(dataCellRegPath, dataCellRegFilename);
    tmp = load(dataCellRegFullFilename);
    CellRegTable = tmp.CellRegTable;
    

    if ~isempty(CellRegTable)
        % Add the registered cell names to each of the cells
        [MapsData] = add_registered_cell_names_to_table(MapsData, CellRegTable);
        [MapsSquareData] = add_registered_cell_names_to_table(MapsSquareData, CellRegTable);
        
        [MapsDataExcluded] = add_registered_cell_names_to_table(MapsDataExcluded, CellRegTable);
        [MapsSquareDataExcluded] = add_registered_cell_names_to_table(MapsSquareDataExcluded, CellRegTable);
    end
    
    % Load the sfp data that has been manually scored
    if analysisSettings.CALCIUM_USE_ONLY_MANUAL_SFP_SCORES
        [dataSfpFilename, dataSfpPath] = uigetfile({'*.mat'}, 'Select calcium_sfp_data with scores mat file.');
        dataSfpFullFilename = fullfile(dataSfpPath, dataSfpFilename);
        tmp2 = load(dataSfpFullFilename);
        SfpData = tmp2.SfpData;
        SfpData(SfpData.sfp_score == 0,:) = [];
        goodCellNames = cell(size(SfpData,1),1);
        for iRow = 1:size(SfpData,1)
           goodCellNames{iRow} = sprintf('%s.t', SfpData.animalSessionCellName{iRow}); 
        end
        % Use only good cells
        MapsData(~ismember(MapsData.cellName, goodCellNames),:) = [];
        MapsSquareData(~ismember(MapsSquareData.cellName, goodCellNames),:) = [];
    end
    
    
    % Store results
    analysisInput = [];
    analysisInput.MapsData = MapsData;
    analysisInput.MapsSquareData = MapsSquareData;
    analysisInput.MapsDataExcluded = MapsDataExcluded;
    analysisInput.MapsSquareDataExcluded = MapsSquareDataExcluded;
    
    %if analysisSettings.CALCIUM_USE_ONLY_REGISTERED_CELLS_IN_N_SESSIONS > 1
        analysisInput.CellRegTable = CellRegTable;
    %end
end % function

