% provides electron density in a Gaussian-shaped DOS as a function of the
% chemical potential and temperature
%
% here, there is no discrimination between donor/acceptor/neutral like
% states

% input: DOS_Gauss_item     .. data set associated to Gaussian DOS
%        chemical_potential .. chemical potential / eV        
%        T                  .. temperature / K
% output: determined density

function [ density ] = GetDensityInLevel(chemical_potential,DOS_level_item,T )
  

    density = 0;

%   occupation = f_FD(EL, chem_pot,T) * N_L
    density = FermiDirac(DOS_level_item.E_ref, chemical_potential,T) * DOS_level_item.N;
end
