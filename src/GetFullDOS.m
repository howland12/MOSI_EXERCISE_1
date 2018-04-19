% this function provides the full density of states in a vector 
% consistent with the energy interval
% input:    DOS_admin .. data containing DOS information
% output:   DOS       .. vector containing DOS

function [DOS] = GetFullDOS(DOS_admin)
    
    DOS = DOS_admin(1).energies;
    max_index = size(DOS_admin,2);

    for i=2:max_index
        DOS = DOS + DOS_admin(i).energies;
    end
end