% provides electron density in a Gaussian-shaped DOS as a function of the
% chemical potential and temperature
%
% here, there is no discrimination between donor/acceptor/neutral like
% states
%
% here, the integration is done independent of the energy mesh

% input: DOS_Gauss_item     .. data set associated to Gaussian DOS
%        chemical_potential .. chemical potential / eV        
%        T                  .. temperature / K
% output: determined density

function [ density ] = GetDensityInGauss(chemical_potential,DOS_Gauss_item,T )

    density = 0;

% check that neither sqrt- or delta-shaped items are considered

    if ((DOS_Gauss_item.label(2) ~= 'L') && ...
        (DOS_Gauss_item.label(2) ~= 'B'))

        func_handle = @(E)FDIntegrantGauss(E,DOS_Gauss_item,...
                                           chemical_potential,T);
        density = integral(func_handle,-inf,inf);

    end
end

