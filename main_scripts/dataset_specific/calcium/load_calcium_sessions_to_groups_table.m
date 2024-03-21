% Marc Normandin, Muzzio Lab, Psychological & Brain Sciences, University of Iowa, 2024
function [SessionsToGroups] = load_calcium_sessions_to_groups_table()
    R(1).animalName = 'CMG129_CA1';
    R(1).sessionName = 's1';
    R(1).dayNum = 1;
    R(1).dayLabel = 'Day 1';
    R(1).groupId = 1;
    R(1).groupLabel = 'Day 1';
    
    R(2).animalName = 'CMG129_CA1';
    R(2).sessionName = 's2';
    R(2).dayNum = 2;
    R(2).dayLabel = 'Day 2';
    R(2).groupId = 2;
    R(2).groupLabel = 'Day 2';
    
    R(3).animalName = 'CMG129_CA1';
    R(3).sessionName = 's3';
    R(3).dayNum = 3;
    R(3).dayLabel = 'Day 3';
    R(3).groupId = 3;
    R(3).groupLabel = 'Day 3';
    
    
    
    
    R(4).animalName = 'CMG154_CA1';
    R(4).sessionName = 's1';
    R(4).dayNum = 1;
    R(4).dayLabel = 'Day 1';
    R(4).groupId = 1;
    R(4).groupLabel = 'Day 1';
    
    R(5).animalName = 'CMG154_CA1';
    R(5).sessionName = 's2';
    R(5).dayNum = 2;
    R(5).dayLabel = 'Day 2';
    R(5).groupId = 2;
    R(5).groupLabel = 'Day 2';
    
    R(6).animalName = 'CMG154_CA1';
    R(6).sessionName = 's3';
    R(6).dayNum = 3;
    R(6).dayLabel = 'Day 3';
    R(6).groupId = 3;
    R(6).groupLabel = 'Day 3';




    R(7).animalName = 'CMG161_CA1';
    R(7).sessionName = 's1';
    R(7).dayNum = 1;
    R(7).dayLabel = 'Day 1';
    R(7).groupId = 1;
    R(7).groupLabel = 'Day 1';
    
    R(8).animalName = 'CMG161_CA1';
    R(8).sessionName = 's2';
    R(8).dayNum = 2;
    R(8).dayLabel = 'Day 2';
    R(8).groupId = 2;
    R(8).groupLabel = 'Day 2';
    
    R(9).animalName = 'CMG161_CA1';
    R(9).sessionName = 's3';
    R(9).dayNum = 3;
    R(9).dayLabel = 'Day 3';
    R(9).groupId = 3;
    R(9).groupLabel = 'Day 3';
    
    
    
    
    
    R(10).animalName = 'CMG162_CA1';
    R(10).sessionName = 's1';
    R(10).dayNum = 1;
    R(10).dayLabel = 'Day 1';
    R(10).groupId = 1;
    R(10).groupLabel = 'Day 1';
    
    R(11).animalName = 'CMG162_CA1';
    R(11).sessionName = 's4';
    R(11).dayNum = 2;
    R(11).dayLabel = 'Day 2';
    R(11).groupId = 2;
    R(11).groupLabel = 'Day 2';
    
    R(12).animalName = 'CMG162_CA1';
    R(12).sessionName = 's5';
    R(12).dayNum = 3;
    R(12).dayLabel = 'Day 3';
    R(12).groupId = 3;
    R(12).groupLabel = 'Day 3';
    
    
    
    
    
    R(13).animalName = 'CMG169_CA1';
    R(13).sessionName = 's1';
    R(13).dayNum = 1;
    R(13).dayLabel = 'Day 1';
    R(13).groupId = 1;
    R(13).groupLabel = 'Day 1';
    
    R(14).animalName = 'CMG169_CA1';
    R(14).sessionName = 's2';
    R(14).dayNum = 2;
    R(14).dayLabel = 'Day 2';
    R(14).groupId = 2;
    R(14).groupLabel = 'Day 2';
    
    R(15).animalName = 'CMG169_CA1';
    R(15).sessionName = 's3';
    R(15).dayNum = 3;
    R(15).dayLabel = 'Day 3';
    R(15).groupId = 3;
    R(15).groupLabel = 'Day 3';

SessionsToGroups = struct2table(R);
    
end % function
