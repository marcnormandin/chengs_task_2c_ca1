function [s] = load_dataset(folder)
    [BestAligned, analysisSettings] = load_best_aligned(folder);
    [MapsData] = load_maps_data(folder, analysisSettings);

    % Verify the integrity of the data because we have to index between the tables.
    [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned);
    if n_errors > 0
        error('Data compromised.\n');
    end

    % Store it in a struct to make it easier to manage
    s = [];
    s.folder = folder;
    s.analysisSettings = analysisSettings;
    s.MapsData = MapsData;
    s.BestAligned = BestAligned;
end

function [n_errors, bad_ba_rows] = verify_data_integrity(MapsData, BestAligned)
    n_errors = 0;
    bad_ba_rows = [];
    for bestAligned_row_index = 1:size(BestAligned,1)
        animalName = BestAligned.animalName{bestAligned_row_index};
        sessionName = BestAligned.sessionName{bestAligned_row_index};
        cellName = BestAligned.cellName{bestAligned_row_index};
        bestAligned_cellDataInds = BestAligned.cellDataInds{bestAligned_row_index};
        [output] = get_cell_data_by_name(MapsData, animalName, sessionName, cellName);
        mapsData_cellDataInds = output.cellDataInds;
        n_inds = length(bestAligned_cellDataInds);
        if n_inds ~= length(intersect(bestAligned_cellDataInds, mapsData_cellDataInds))
            %fprintf('Error with best aligned row index: %d\n', bestAligned_row_index);
            n_errors = n_errors + 1;
            bad_ba_rows(n_errors) = bestAligned_row_index;
        end
    end % bestAligned_row_index

    if ~isempty(bad_ba_rows)
        error('Corrupt data')
    end
end


function [MapsData] = load_maps_data(folder, analysisSettings)
    tmp = load(fullfile(folder, 'analysis_input.mat'));
    MapsData = tmp.analysisInput.MapsSquareData;
    eliminateUnmatched = false;
    MapsData = helper_add_group_labels_to_table(MapsData, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function


function [BestAligned, analysisSettings] = load_best_aligned(folder)
    tmp = load(fullfile(folder, 'analysis_results.mat'));
    analysisSettings = tmp.analysisSettings;
    BestAligned = tmp.analysisResults.BestAligned;
    eliminateUnmatched = false;
    BestAligned = helper_add_group_labels_to_table(BestAligned, load_sessions_to_groups_table(analysisSettings), eliminateUnmatched);
end % function
