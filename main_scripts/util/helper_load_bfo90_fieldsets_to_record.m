% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [fieldSets] = helper_load_bfo90_fieldsets_to_record()
    fieldSets(1).name = 'all';
    fieldSets(1).fieldNames = {'any'};

    fieldSets(2).name = 'context1';
    fieldSets(2).fieldNames = {'context1'};

    fieldSets(3).name = 'context2';
    fieldSets(3).fieldNames = {'context2'};

    fieldSets(4).name = 'within';
    fieldSets(4).fieldNames = {'context1', 'context2'};

    fieldSets(5).name = 'across';
    fieldSets(5).fieldNames = {'different'};
end
