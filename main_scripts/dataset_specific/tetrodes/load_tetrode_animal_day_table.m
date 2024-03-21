% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [R] = load_tetrode_animal_day_table()
    R = [];

    day1.animalNames = {'K1_CA1', 'AK42_CA1', 'AK74_CA1', 'CMG159_recut', 'HGY1_recut', 'JJ9_CA1', 'MG1'};
    day1.sessionNames = {'d1',       'd7',         'd1',         's6',         'd1',      'd1',      'd6'};

    day2.animalNames = {'K1_CA1', 'AK42_CA1', 'AK74_CA1', 'CMG159_recut',               'JJ9_CA1', 'MG1'};
    day2.sessionNames = {'d2',       'd8',        'd2',          's7',                   'd2',       'd7'};

    day3.animalNames = {'K1_CA1', 'AK42_CA1', 'AK74_CA1', 'CMG159_recut',               'JJ9_CA1', 'MG1', };
    day3.sessionNames = {'d4',       'd9',         'd3',          's8',                 'd3',     'd11',   };


    % process a day
    %days = [];
    clear days
    days(1) = day1;
    days(2) = day2;
    days(3) = day3;

    k = 1;
    for iDay = 1:length(days)
        day = days(iDay);
        for iRec = 1:length(day.animalNames)
            R(k).animalName = day.animalNames{iRec};
            R(k).sessionName = day.sessionNames{iRec};
            R(k).dayNum = iDay;
            R(k).dayLabel = sprintf('Day %d', R(k).dayNum);
            k = k + 1;
        end
    end

    R = struct2table(R);
end % function


