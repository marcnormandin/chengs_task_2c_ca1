% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [hFig] = plot_bestaligned_registered(BestAligned, MapsData, registeredCellName, doNormalization)
    hFig = figure('position', get(0, 'screensize'));

    maxTrials = 6+2;


    cdata = BestAligned(ismember(BestAligned.registeredCellName, registeredCellName),:);
    numSessions = size(cdata,1);

    %mm = reshape(1:(maxTrials*numSessions*2),2*numSessions,maxTrials)';

    for iSession = 1:numSessions
        cellName = cdata.cellName{iSession};

        animalName = cdata.animalName{iSession};
        sessionName = cdata.sessionName{iSession};
        [output] = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        maps1 = output.maps(:,:,output.contextIds == 1);
        maps2 = output.maps(:,:,output.contextIds == 2);
        rot1 = cdata.rotationSequence1{iSession};
        rot2 = cdata.rotationSequence2{iSession};
        
        % Rotate the maps to what we found was best
        for i1 = 1:length(rot1)
           m1 = maps1(:,:,i1);
           if rot1(i1) == 1
               m1 = rot90(m1,2);
           end
           
           if doNormalization
               m1 = m1 ./ sum(m1, 'all');
           end
               
           maps1(:,:,i1) = m1;
           

        end
        for i2 = 1:length(rot2)
           m2 = maps2(:,:,i2);
           if rot2(i2) == 1
               m2 = rot90(m2,2);
           end
           
           if doNormalization
               m2 = m2 ./ sum(m2, 'all');
           end
           
           maps2(:,:,i2) = m2;
        end
        
        mean1 = mean(maps1, 3);
        error1 = std(maps1, 0, 3) ./ sqrt(size(maps1,3));
        
        meanMean1 = mean(mean1, 'all');
        meanError1 = mean(error1, 'all');
        
        mean2 = mean(maps2, 3);
        error2 = std(maps2, 0, 3) ./ sqrt(size(maps2,3));
        
        meanMean2 = mean(mean2, 'all');
        meanError2 = mean(error2, 'all');
        

        PP = 2*numSessions;
        QQ = maxTrials;

        
        
        for i1 = 1:length(rot1)
            k1 = 2*maxTrials*(iSession-1) + i1;
           subplot(PP,QQ,k1);
           m1 = maps1(:,:,i1);
%            if rot1(i1) == 1
%                m1 = rot90(m1,2);
%            end
           ml_imagesci( m1 );
           axis equal tight off
           title(sprintf('%s Context 1', sessionName));
           colormap jet
           colorbar
        end
        for i2 = 1:length(rot2)
            k2 = 2*maxTrials*(iSession-1) + maxTrials + i2;
           subplot(PP,QQ,k2);
           m2 = maps2(:,:,i2);
%            if rot2(i2) == 1
%                m2 = rot90(m2,2);
%            end
           ml_imagesci( m2 );
           axis equal tight off
           title(sprintf('%s Context 2', sessionName)); 
           colormap jet
           colorbar
        end

        % Show the best aligned maps
        subplot(PP,QQ,2*maxTrials*(iSession-1) + maxTrials-1);
        %ml_imagesci( cdata.averageMapContext1{iSession} );
        ml_imagesci( mean1 );
        axis equal tight off
        title(sprintf('%s Mean Map\nContext 1\nMean pixel mean %0.2e', sessionName, meanMean1));
        colormap jet
        %colorbar
        
        subplot(PP,QQ,2*maxTrials*(iSession-1) + maxTrials);
        %ml_imagesci( cdata.averageMapContext1{iSession} );
        ml_imagesci( error1 );
        axis equal tight off
        title(sprintf('%s Error Map\nContext 1\nMean pixel error %0.2e', sessionName, meanError1));
        colormap jet
        %colorbar

        subplot(PP,QQ,2*maxTrials*(iSession-1) + 2*maxTrials-1);
        %ml_imagesci( cdata.averageMapContext2{iSession} );
        ml_imagesci( mean2 );
        axis equal tight off
        title(sprintf('%s Mean Map\nContext 2\nMean pixel mean %0.2e', sessionName, meanMean2));
        colormap jet
        %colorbar
        
        subplot(PP,QQ,2*maxTrials*(iSession-1) + 2*maxTrials);
        %ml_imagesci( cdata.averageMapContext2{iSession} );
        ml_imagesci( error2 );
        axis equal tight off
        title(sprintf('%s Error Map\nContext 2\nMean pixel error %0.2e', sessionName, meanError2));
        colormap jet
        %colorbar
    end
    if doNormalization
        sgtitle(sprintf('Registred Cell (NORMALIZED MAPS)\n%s', registeredCellName), 'interpreter', 'none')
    else
        sgtitle(sprintf('Registred Cell (NO NORMALIZATION)\n%s', registeredCellName), 'interpreter', 'none')
    end
end % function
