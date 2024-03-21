% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [M] = get_match_table_by_sessions(TableA, TableB)
    % Get matching sessions between the two tables

    matchingRows = ml_util_find_row_matches([TableA.animalName, TableA.sessionName], [TableB.animalName, TableB.sessionName]);
    M = [];
    for iMatch = 1:size(matchingRows,1)
        M(iMatch).animalName = TableA.animalName{matchingRows(iMatch,1)}; % should be the same as for B
        M(iMatch).sessionName = TableA.sessionName{matchingRows(iMatch,1)};
       
        M(iMatch).TableARowIndex = matchingRows(iMatch,1);
        M(iMatch).TableBRowIndex = matchingRows(iMatch,2);
    end
    if ~isempty(M)
        M = struct2table(M);
    end
end
