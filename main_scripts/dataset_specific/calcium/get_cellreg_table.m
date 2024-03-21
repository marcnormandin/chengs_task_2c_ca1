% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [C] = get_cellreg_table(cellRegFolder)
    % We need to load the cellreg results, and map them back into sessions
    % tables
    [cellToIndexMap, cellScores] = load_cellreg_results(cellRegFolder);

    files = dir(fullfile(cellRegFolder, 'sfp_*.mat'));


    % Be careful about sorting because if we have more than 9 this will mess up
    if length(files) > 9
        error('Files will need to be properly sorted.\n');
    end

    T = [];
    for iFile = 1:length(files)
        f = files(iFile);

        fn = fullfile(f.folder, f.name);
        data = load(fn);
        T(iFile).animalName = data.animalName;
        T(iFile).sessionName = data.sessionName;
        T(iFile).SFP = data.SFP;
    end
    T = struct2table(T);
    % The SFPs are #Cells x Video Height x Video Width arrays

    numSessions = size(cellToIndexMap,2);
    if numSessions ~= length(files)
        error('Number of sessions do not match!\n');
    end

    C = [];
    k = 1;
    for iGlobal = 1:size(cellToIndexMap,1)
         globalId = iGlobal;
         for iSession = 1:numSessions
            cid = cellToIndexMap(iGlobal, iSession);
            cellScore = cellScores(iGlobal);
            if cid ~= 0
                C(k).animalName = T.animalName{iSession};
                C(k).sessionName = T.sessionName{iSession};
                C(k).sessionCellId = cid;
                C(k).cellName = sprintf('%s_%s_%d.t', C(k).animalName, C(k).sessionName, cid);
                
                sfps = T.SFP{iSession};
                C(k).spatialFootprint = squeeze( sfps(cid,:,:) );
                
                C(k).globalCellId = globalId;
                C(k).registeredCellName = sprintf('%s_%d', C(k).animalName, C(k).globalCellId);
                C(k).cellScore = cellScore;
                
                k = k + 1;
            end
         end
    end
    C = struct2table(C);


    % Counts
    globalCellIdCounts = histcounts(C.globalCellId, 1:max(C.globalCellId)+1);
    globalCellIdCountsByRow = zeros(size(C,1),1);
    for i = 1:size(C,1)
        globalCellIdCountsByRow(i) = globalCellIdCounts(C.globalCellId(i));
    end

    C.sharedCount = globalCellIdCountsByRow;
end % function 

function [cellToIndexMap, cellScores] = load_cellreg_results(cellRegFolder)
    files = dir(fullfile(cellRegFolder, 'cellRegistered*.mat'));
    if length(files) ~= 1
        error('No file found.\n');
    end
    
    fn = fullfile(files(1).folder, files(1).name);
    tmp = load(fn);
    
    cellToIndexMap = tmp.cell_registered_struct.cell_to_index_map;
    cellScores = tmp.cell_registered_struct.cell_scores;
end
