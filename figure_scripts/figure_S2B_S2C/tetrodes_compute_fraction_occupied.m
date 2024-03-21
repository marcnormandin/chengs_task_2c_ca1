% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [C] = tetrodes_compute_fraction_occupied(ratemaps, useSmoothed)
    C = [];
    k = 1;
    numAnimals = length(ratemaps.data);
    for iAnimal = 1:numAnimals
        d1 = ratemaps.data(iAnimal);
        animalName = d1.animal;
        groupName = d1.group;
        numSessions = length(d1.session);
        for iSession = 1:numSessions
            d2 = d1.session(iSession);
            sessionName = d2.name;
    
            numTrials = length(d2.trial);
            for iTrial = 1:numTrials
                d3 = d2.trial(iTrial);
                trialId = d3.trialNum;
                if useSmoothed
                    occupancyMap = d3.smoothedTime > 0.05;
                else
                    occupancyMap = d3.timeMap > 0.05;
                end
                
    
                % compute
                % no padding cuz no miniscope
                R = compute_whatever_we_need_per_map(occupancyMap, 0);
            
                C(k).animalName = animalName;
                C(k).groupName = groupName;
                C(k).sessionName = sessionName;
                C(k).trialId = trialId;
                C(k).occupancyMap = occupancyMap;

                fieldNames = fieldnames(R);
                for iField = 1:length(fieldNames)
                    fieldName = fieldNames{iField};
                    C(k).(fieldName) = R.(fieldName);
                end
                k = k + 1;
            end
        end
    end
    C = struct2table(C);
end % function

