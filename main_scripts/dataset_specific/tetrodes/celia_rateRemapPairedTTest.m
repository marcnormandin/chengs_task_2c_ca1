% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024

%% Paired TTest using cell as statistical unit.

%load('../output/ratemaps-2021_10_1_12_35_44.mat');

MIN_TRIALS = 4;

% animalName = 'AK42_CA1';
% sessionName = 'd7';
% cellName = 'TT2_02.t';

% animalName = 'CMG159_recut';
% sessionName = 's3';
% cellName = 'CMG159-s3-TT2_4.t';

animalName = 'AK74_CA1';
sessionName = 'd2';
cellName = 'TT4_06.t';


[withinMean, acrossMean, results, isValid] =  ratemaps_within_across_mean_for_cell_compute(MIN_TRIALS, ratemaps, animalName, sessionName, cellName)





%% per day paired ttest
avg = nan(3,2);
sem = nan(3,2);
for d = 1:3
    withinDay = cell2mat(within(:,d));
    acrossDay = cell2mat(across(:,d));
    
    %[H,P] = ttest(withinDay, acrossDay, 'tail', 'left');
    [H,P] = ttest(withinDay, acrossDay);
    avg(d,:) = [nanmean(withinDay), nanmean(acrossDay)];
    sem(d,:) = [nanstd(withinDay) ./ sqrt(sum(~isnan(withinDay))), nanstd(acrossDay) ./ sqrt(sum(~isnan(acrossDay)))];
    fprintf('Two Tailed paired ttest for rate remapping day %d: \n \t p = %.3f\n', d,P);
end

% Example data
b = bar(avg, 'grouped');
hold on
% Find the number of groups and the number of bars in each group
[ngroups, nbars] = size(avg);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, avg(:,i), sem(:,i), 'k', 'linestyle', 'none');
end
hold off


function [withinMean, acrossMean, results, isValid] =  ratemaps_within_across_mean_for_cell_compute(MIN_TRIALS, ratemaps, animalName, sessionName, cellName)
    animalNames = unique(extractfield(ratemaps.data, 'animal'));
    indAnimal = find(ismember(animalNames, animalName));
    if isempty(indAnimal)
        error('Animal (%s) not found\n', animalName);
    end

    % session names available (NOT DAYS)
    sessionNames = {ratemaps.data(indAnimal).session.name};
    sessionInd = find(ismember(sessionNames, sessionName));
    if isempty(sessionInd)
        error('Session %s not found for %s\n', sessionName, animalName);
    end

    sessionData = ratemaps.data(indAnimal).session(sessionInd);
    trialData = ratemaps.data(indAnimal).session(sessionInd).trial;

    %cellNames
    cellNames = {};
    numTrials = length(trialData);
    for iTrial = 1:numTrials
        labels = trialData(iTrial).label;
        for i = 1:length(labels)
            cellNames{length(cellNames)+1} = labels{i};
        end
    end
    cellNames = unique(cellNames);

    cellInd = find(ismember(cellNames, cellName));
    if isempty(cellInd)
        error('Cell %s not found\n', cellName)
    end
    trials = extractfield(trialData, 'trialNum');
    contexts = extractfield(trialData, 'context');
    digs = extractfield(trialData, 'dig');


    withinMean = [];
    acrossMean = [];
    isValid = false;
    results = [];

    mfr = calculateMfrVector(sessionData, cellName);

    trialsAvail = trials(ismember(digs, {'Corr', 'Geo', 'Feat', 'Wrong'}) & ~isnan(mfr));
    if length(trialsAvail) < MIN_TRIALS
        return;
    end
    
    isValid = length(unique(contexts(trialsAvail))) >= 2;
    if ~isValid
        return;
    end

    com = nchoosek(trialsAvail,2);
    numComparisons = size(com,1);
    totalWithin = nan(numComparisons,1);
    totalAcross = nan(numComparisons,1);
    mfrdiff = nan(numComparisons,1);
    for iComparison = 1:numComparisons

        trialA = com(iComparison,1);
        trialB = com(iComparison,2);

        mfrdiff(iComparison) = abs(mfr(trialA) - mfr(trialB));

        if contexts(trialA) == contexts(trialB)
            totalWithin(iComparison) = mfrdiff(iComparison);
        else
            totalAcross(iComparison) = mfrdiff(iComparison);
        end
    end
    
    results = [com, mfrdiff];

    withinMean = nanmean(totalWithin);
    acrossMean = nanmean(totalAcross);

end % function
