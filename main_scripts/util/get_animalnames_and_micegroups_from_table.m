% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = get_animalnames_and_micegroups_from_table(T)
% This function returns a table that contains the unique pairs of
% animalName and mouseGroup in the table T.

uniqueAnimalNames = unique(T.animalName);
uniqueMouseGroups = unique(T.mouseGroup);

nRows = size(T,1);
H = nan(nRows,2);

for i = 1:length(uniqueAnimalNames)
    animalName = uniqueAnimalNames{i};
    H(ismember(T.animalName, animalName),1) = i;
end

for i = 1:length(uniqueMouseGroups)
    mouseGroup = uniqueMouseGroups{i};
    H(ismember(T.mouseGroup, mouseGroup),2) = i;
end
u = unique(H, 'rows');

h = [];
for i = 1:size(u,1)
   h(i).animalName = uniqueAnimalNames{u(i,1)};
   h(i).mouseGroup = uniqueMouseGroups{u(i,2)};
end
R = struct2table(h);
end % function
