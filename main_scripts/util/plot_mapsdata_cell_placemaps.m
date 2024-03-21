% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024

% Set these
%MapsData = analysisInput.MapsData;
%MapsDataInfo = get_sessions_and_cellnames_from_table(MapsData);


%%
close all


% Select a cell
iRow = 3004;
%inds = find(ismember(MapsDataInfo.animalSessionCellName, {'CMG169_CA1_s3_CMG169_CA1_s3_219.t'}));
%inds = find(ismember(MapsDataInfo.cellName, {'CMG169_CA1_s3_219.t'}));
%iRow = inds(1);

animalName = MapsDataInfo.animalName{iRow};
sessionName = MapsDataInfo.sessionName{iRow};
cellName = MapsDataInfo.cellName{iRow};

% Plot the cell
plot_cells_maps(MapsData, animalName, sessionName, cellName)

function plot_cells_maps(MapsData, animalName, sessionName, cellName)
    MAXIMUM_TOTAL_TRIALS = 12;

    % Get the cell data
    cdata = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);

    [subplotTrialIds, ~] = helper_subplot_indices(MAXIMUM_TOTAL_TRIALS, cdata.trialIds, cdata.contextIds);

    % Plot the cell's maps
    nMaps = cdata.numEntries;
    for iMap = 1:nMaps
        tid = cdata.trialIds(iMap);
        cid = cdata.contextIds(iMap);

        sid = find(subplotTrialIds == tid);

        subplot(2,MAXIMUM_TOTAL_TRIALS/2, sid)
        m = squeeze(cdata.maps(:,:,iMap));
        imagesc(m)
        axis equal tight off
        set(gca, 'ydir', 'reverse')
        title(sprintf('T%d C%d', tid, cid))
        colormap jet
    end % iMap
    sgtitle(sprintf('%s %s %s', animalName, sessionName, cellName), 'interpreter', 'none')
end % function
    

function [rateMatrixTrialIds, rateMatrixContextIds] = helper_subplot_indices(rateMatrixNumTrials, trialIds, contextIds)
%     rateMatrixContextIds = repmat([1,2], rateMatrixNumTrials/2,1); % will only work if context is 1 or 2
%     rateMatrixContextIds = rateMatrixContextIds(:); % 1 1 1 1 1 1 2 2 2 2 2 2
%     rateMatrixContextIds = reshape(rateMatrixContextIds, 1, length(rateMatrixContextIds));
    rateMatrixContextIds = [ones(1,rateMatrixNumTrials/2), 2*ones(1, rateMatrixNumTrials/2)];

    
    % Sort the rate matrix rows and cols by context, so that we can see
    % patterns based on context more easily.
    % Get the unique combinations of (trial id, context id)
    urows = unique([trialIds, contextIds], 'rows');
    [surows, ~] = sortrows(urows, 2);
    
    context1Pairs = surows(surows(:,2)==1,:);
    context2Pairs = surows(surows(:,2)==2,:);
    
    hasContext1 = ~isempty(context1Pairs);
    hasContext2 = ~isempty(context2Pairs);
    
    if hasContext1
        % Make sure that all trial ids are uniformly even or odd
        tids = context1Pairs(:,1);
        context1IsEven = all(mod(tids,2)==0);
    end
    
    if hasContext2
        % Make sure that all trial ids are uniformly even or odd
        tids = context2Pairs(:,1);
        context2IsEven = all(mod(tids,2)==0);
    end
    
    % possibililty A: Context 1 has odd-valued trials and Context 2 has
    % even-valued trials.
    rateMatrixTrialIdsA = [1:2:rateMatrixNumTrials, 2:2:rateMatrixNumTrials];
        
    % possibility B: Context 1 has even-valued trials and Context 2 has
    % odd-valued trials.
    rateMatrixTrialIdsB = [2:2:rateMatrixNumTrials, 1:2:rateMatrixNumTrials];

    rateMatrixTrialIds = nan;
    
    if hasContext1 && hasContext2
        if context1IsEven && ~context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        elseif ~context1IsEven && context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        else
            error('Invalid contexts do not contain uniformly even or odd trial ids')
        end
    elseif hasContext1
        if context1IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        else
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        end
    elseif hasContext2
        if context2IsEven
            rateMatrixTrialIds = rateMatrixTrialIdsA;
        else
            rateMatrixTrialIds = rateMatrixTrialIdsB;
        end
    else
        error('Something is messed up with the context ids')
    end
end % function
