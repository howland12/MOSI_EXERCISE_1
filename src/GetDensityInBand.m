% Here the charge density within the conduction or valence band is 
% computed with the approximated Fermi integral 
%
%
function [charge_density] = ...
    GetDensityInBand(chemical_potential,E_edge, eff_mass,T)

    k_eV = 1.3806504E-23 / 1.602176487E-19; % Boltzmann constant in eV/K
    arg = -abs(chemical_potential-E_edge)/k_eV/T;

    charge_density = 2/ sqrt(pi) * DensityOfBandStates(eff_mass,T) * ...
                    FIntegrationForBands(arg);

end