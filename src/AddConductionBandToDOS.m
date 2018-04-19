% adds DOS of conduction band to a given DOS
% provides vector to store the hole occupation

% input: DOS .. previous DOS vector
%        E   .. energy vector
%        E_C .. minimum of conduction band   (eV)
%        eff_mass .. effective electron mass (in units of me)

% output: DOS vector with number of states per energy interval 
%         corresponding to requested Gaussian added
%         occupation vector (zero-valued)


function [ DOS_admin] = AddConductionBandToDOS( DOS_admin, E, E_C, eff_mass)


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

energies_above_E_C = find(E>E_C);



    DOS = zeros(size(E),'like',E);

    DOS(energies_above_E_C) = 1 * prefactor * (eff_mass)^(3/2).* sqrt(E(energies_above_E_C)-E_C);
   
    DOS_admin = AddContribtionToDOS(DOS_admin, DOS, 'CB','N',E_C,0,0);

%occ_vector = zeros(size(E),'like',E);
end



