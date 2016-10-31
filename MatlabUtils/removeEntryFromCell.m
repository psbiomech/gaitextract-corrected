% function newCell = removeEntryFromCell(oldCell, entry)
% The entry must be a string.
% This function removes all instances of the entry
% (even if it appears more than once)
% Tim Dorn
% September 2010
% ----------------------------------------------------------

function newCell = removeEntryFromCell(oldCell, entry)

newCell = {};
l = length(oldCell);

for i = 1:l
    if ~strcmp(oldCell{i}, entry)
        newCell = [newCell, oldCell{i}];
    end
end
