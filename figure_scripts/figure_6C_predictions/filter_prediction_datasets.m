% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [filteredPredictionDatasets] = filter_prediction_datasets(predictionDatasets, params)
    filteredPredictionDatasets = cell(length(predictionDatasets),1);
    for iDataset = 1:length(predictionDatasets)
        predictionDataset = predictionDatasets{iDataset};
        filteredPredictionDataset = filter_prediction_dataset(predictionDataset, params);
        filteredPredictionDatasets{iDataset} = filteredPredictionDataset;
    end % iDataset
end % function

