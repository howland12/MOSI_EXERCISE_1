% calculate the effective density of states for a given effective mass and
% temperature

% input:    eff_mass  .. effective mass (units of electron mass me)
%           T         .. temperature (K)
% output:   N         .. effective density of states for given temperature
%                     ..   in cubic meter (!)

function [N] = DensityOfBandStates(eff_mass,T)

%   considering the parameters
%
%   me = 9.11e-31; % kg
%   k = 1.38e-23;  % SI
%   h = 6.626E-34; % SI
%
%   one arrives at a prefactor of
%
%   prefactor = 2*(2*pi*k*me)^(3/2)/h^3;

    prefactor = 4.8266e+21;

    N = prefactor * (eff_mass * T)^(3/2);

end