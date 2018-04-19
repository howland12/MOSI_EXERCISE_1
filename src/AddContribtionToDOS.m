% stores information on new DOS contribution into DOS_admin
%
% input  DOS_admin .. previous DOS_admin
%        DOS       .. vector containing the same number of bins as energy
%                     interval
%        label     .. 'D', 'A','em' 
%                     will be assigned by AddXXXToDOS 
%                               (XXX <> 'contribution')
%        type      .. 'N', 'P' 
%                     indicator whether states become negatively or  
%                     positively charged upon filling
%        E_ref     .. position of contribution in energy interval (in eV)
%        param     .. optional additional parameter, e.g., Gaussian width
%        N_DOS     .. total density of states (in inverse cubic meters)
%       

% output DOS_admin .. updated DOS_admin

function [DOS_admin] = AddContribtionToDOS(DOS_admin, DOS, label,type,...
                                           E_ref,param,N_DOS)

    % check for valid entries in label and type

    tempDOS.energies = DOS;
    tempDOS.label = label;
    tempDOS.type = type;

    tempDOS.E_ref = E_ref;
    tempDOS.param = param;
    tempDOS.N = N_DOS;

    % is this the first entry in DOS?
    % yes: overwrite dummy entry
    % no: expand list and set values

    index = size(DOS_admin,2);

    % if ((index == 1)&&(DOS_admin(1).label ~= 'empty') )
    if (DOS_admin(index).label == 'em')
        DOS_admin(index) = tempDOS;        
    else
        DOS_admin = [DOS_admin,tempDOS];
    end   

end