% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
% This loads celias tetrode data into a form that can be used by the
% population vector code
function [T] = load_celia_data_structure_as_table_v2(data)
    %data = load(dataFilename);
    %data = data.ratemaps;

    animalNames = extractfield(data, 'animal');
    numAnimals = length(animalNames);
    %animalIds = 1:numAnimals;

    
    T = [];
    k = 1;
    for iAnimal = 1:numAnimals
        sessionsData = data(iAnimal).session;
        sessionNames = extractfield(sessionsData, 'name');
        animalName = data(iAnimal).animal;

        for iSession = 1:length(sessionNames)
            sessionName = sessionNames{iSession};
            trialData = sessionsData(iSession).trial;

            numTrials = length(trialData);

            % Get all of the unique cell names that are present for the
            % SESSION so we can give each a unique integer id
            cellNames = [];
            for iTrial = 1:numTrials
                tdata = trialData(iTrial);
                for iCell = 1:length(tdata.label)
                   %cellNames = cat(1, cellNames, tdata.label{iCell});
                   kk = length(cellNames)+1;
                   cellNames{kk} = tdata.label{iCell};
                end
            end
            uniqueCellNames = unique(cellNames);
            cellInds = arrayfun(@(x)(find(ismember(cellNames, uniqueCellNames{x}))), 1:length(uniqueCellNames), 'UniformOutput', false);
            uniqueCellIds = nan(length(cellNames), 1);
            for iCell = 1:length(cellInds)
                uniqueCellIds( cellInds{iCell} ) = iCell;
            end
            
            cnk = 1;
            for iTrial = 1:numTrials
                tdata = trialData(iTrial);

                
                
                tid = tdata.trialNum;
                cid = tdata.context;
                numCells = length(tdata.label);
                for iCell = 1:numCells
                   T(k).animalName = animalName;
                   T(k).sessionName = sessionName;
                   T(k).dayNum = iSession;
                   T(k).cellName = tdata.label{iCell};
                   T(k).cellId = uniqueCellIds(cnk);
                   T(k).trialId = tid;
                   T(k).contextId = cid;
                   T(k).map = tdata.maps{iCell};
                   
                   T(k).dig = tdata.dig;
                   
                   T(k).animalSessionCellName = sprintf('%s_%s_%s', T(k).animalName, T(k).sessionName, T(k).cellName);
                                      
                   % The square maps structure doesn't have all of the same
                   % field as the normal ratemaps structure for rectangular
                   % maps so we have to check.
                   if isfield(tdata, 'mfr')
                    T(k).mfr = tdata.mfr{iCell};
                   end
                   
                   if isfield(tdata, 'pfr')
                    T(k).pfr = tdata.pfr{iCell};
                   end
                   
                   if isfield(tdata, 'infoRate')
                    T(k).infoRate = tdata.infoRate{iCell};
                   end
                   
                   if isfield(tdata, 'infoSpike')
                    T(k).infoSpike = tdata.infoSpike{iCell};
                   end
                   
                   k = k + 1;
                   cnk = cnk + 1;
                end
            end
        end
    end
    T = struct2table(T);
end % function

