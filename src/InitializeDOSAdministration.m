 function [DOS_admin] = InitializeDOSAdministration (energy_interval);

    clear DOS_admin;
    DOS_admin.energies = energy_interval;
    DOS_admin.label = 'em'; % em = empty
    DOS_admin.type = 'N';

    DOS_admin.E_ref = 0 ;
    DOS_admin.param = 0;
    DOS_admin.N = 0;  % total density or effective density of states
end