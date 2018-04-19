% adds DOS of valence band to a given DOS
% provides vector to store the hole occupation

% input: E   .. energy vector
%        E_V .. maximum of valence band   (eV)
%        eff_mass .. effective hole mass (in units of me)

% output: entry to DOS with
%         vector with number of states per energy interval 
%         corresponding to requested Gaussian 
%         labeled as VB = valence band
%         type = 'P' for positive charges


function [ DOS_admin ] = AddValenceBandToDOS(DOS_admin,E, E_V, eff_mass)


% considering the parameters
   me = 9.11e-31; % kg
   h = 6.626E-34; % SI
   q = 1.602176565e-19; % SI needed to convert energies from eV to Si
% one arrives at a prefactor

%   when energies given in SI units (J)
%      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)  ;
%   when energies given in eV
%      prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2)*q^(3/2) ;
%
%   this is, because we would like to provide energies in eV and
%   plot the resulting DOS with respect to energy intervals given in eV
%
%   dZ/dE |_eV = dZ/dE |_SI * (dE|_SI)/(dE|_eV) = q dZ/dE |_SI
%   sqrt(E|_SI) = sqrt ( q E|_eV) = sqrt(q) * sqrt(E|_eV)


    prefactor = 8*pi*sqrt(2)*h^(-3)*me^(3/2) * q^(3/2) ;  % bart, pierret,hadley
    %prefactor = 4.2524e+46;

    energies_below_E_V = find(E < E_V);

    DOS = zeros(size(E),'like',E);
    DOS(energies_below_E_V) = 1 * prefactor * (eff_mass)^(3/2).* sqrt(E_V - E(energies_below_E_V));

    DOS_admin = AddContribtionToDOS(DOS_admin, DOS, 'VB','P',E_V,0,0);

    %occ_vector = zeros(size(E),'like',E);
end



