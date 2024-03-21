% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [T] = get_bestaligned_registered_rotated_maps(BestAligned, MapsData, registeredCellName, doNormalization)

    cdata = BestAligned(ismember(BestAligned.registeredCellName, registeredCellName),:);
    numSessions = size(cdata,1);

    %mm = reshape(1:(maxTrials*numSessions*2),2*numSessions,maxTrials)';
    
    T = [];

    for iSession = 1:numSessions
        cellName = cdata.cellName{iSession};

        animalName = cdata.animalName{iSession};
        sessionName = cdata.sessionName{iSession};
        [output] = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        maps1 = output.maps(:,:,output.contextIds == 1);
        maps2 = output.maps(:,:,output.contextIds == 2);
        rot1 = cdata.rotationSequence1{iSession};
        rot1
        rot2 = cdata.rotationSequence2{iSession};
        rot2
        
        disp(length(rot1))
        disp(length(rot2))
        
        
        % Rotate the maps to what we found was best
        for i1 = 1:length(rot1)
           m1 = maps1(:,:,i1);
           if rot1(i1) == 1
               m1 = rot90(m1,2);
           end
           
           if doNormalization
               m1 = m1 ./ nansum(m1, 'all');
           end
               
           maps1(:,:,i1) = m1;
           

        end
        for i2 = 1:length(rot2)
           m2 = maps2(:,:,i2);
           if rot2(i2) == 1
               m2 = rot90(m2,2);
           end
           
           if doNormalization
               m2 = m2 ./ nansum(m2, 'all');
           end
           
           maps2(:,:,i2) = m2;
        end
        
        mean1 = nanmean(maps1, 3);
        error1 = nanstd(maps1, 0, 3) ./ sqrt(size(maps1,3));
        
        meanMean1 = nanmean(mean1, 'all');
        meanError1 = nanmean(error1, 'all');
        
        mean2 = nanmean(maps2, 3);
        error2 = nanstd(maps2, 0, 3) ./ sqrt(size(maps2,3));
        
        meanMean2 = nanmean(mean2, 'all');
        meanError2 = nanmean(error2, 'all');
        
        k = iSession;
        T(k,1).animalName = animalName;
        T(k,1).sessionName = sessionName;
        T(k,1).registeredCellName = registeredCellName;
        T(k,1).sessionCellName = cellName;
        T(k,1).trialIds1 = output.trialIds(output.contextIds == 1);
        T(k,1).trialIds2 = output.trialIds(output.contextIds == 2);
        T(k,1).rotationSequence1 = cdata.rotationSequence1(iSession);
        T(k,1).rotationSequence2 = cdata.rotationSequence2(iSession);
        
        T(k,1).bestAlignedMaps1 = maps1;
        T(k,1).bestAlignedMaps2 = maps2;
        
        T(k,1).meanMap1 = mean1;
        T(k,1).errorMap1 = error1;
        T(k,1).meanMeanMap1 = meanMean1;
        T(k,1).meanErrorMap1 = meanError1;
        
        T(k,1).meanMap2 = mean2;
        T(k,1).errorMap2 = error2;
        T(k,1).meanMeanMap2 = meanMean2;
        T(k,1).meanErrorMap2 = meanError2;
    end % iSession
    
    fprintf('Num sessions = %d\n', numSessions);
    
    %T = struct2table(T);
    
    if numSessions == 1
        T = [];
        return;
    else
        T = struct2table(T);
    end
    %T = struct2table(T);
end % function
