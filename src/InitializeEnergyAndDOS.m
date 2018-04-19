function[EnergyInt,DOS, occup_vec] = InitializeEnergyAndDOS(E_min, E_max,...
                                     E_V, E_C, resolution)

% provides column vectors for 
%        energies
%        DOS (returned zero-valued)
% input: 
%        Emin .. lowest considered energy
%        Emax .. highest considered energy
%        E_V  .. energy of valence band maximum
%        E_C  .. energy of conduction band minimum
%        resolution .. number of energy intervals

% equally spaced energy interval
% EnergyInt = linspace(E_min, E_max, resolution).';

% produce an EnergyInterval being logarithmically divided
number_of_gap_points = round(resolution*(E_C-E_V)/(E_max-E_min));
number_of_val_points = round(resolution*(E_V-E_min)/(E_max-E_min));
number_of_con_points = round(resolution*(E_max-E_C)/(E_max-E_min));

EnergyInt_V = E_V - (E_V-E_min)*logspace(-1,-5,number_of_val_points);
EnergyInt_Gap = linspace(E_V, E_C, number_of_gap_points);
EnergyInt_C = E_C + (E_max-E_C)*logspace(-5,-1,number_of_con_points);

EnergyInt = [EnergyInt_V,EnergyInt_Gap,EnergyInt_C].';

DOS = zeros(size(EnergyInt),'like',EnergyInt);
occup_vec = zeros(size(EnergyInt),'like',EnergyInt);

end

