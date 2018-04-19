% provides the integral kernel of the Fermi Dirac Integral with respect to
% a Gaussian-shaped DOS
%
% note that also 
%         - non-delta (i.e., sharp levels) or 
%         - sqrt-shaped bands 
% can be dealt with this evaluation
%
% input :  E                  .. energy / eV
%          DOS_Gauss_item     .. one data set from DOS admin describing either
%                                a Gauss- or differently shaped DOS (excluding
%                                delta and sqrt-shaped contributions)
%          chemical_potential .. chemical potential / eV
%          T                  .. temperature / K
%   

function [ int_kernel ] = FDIntegrantGauss( E, DOS_Gauss_item, chemical_potential,T)

    int_kernel = 0;
    int_kernel = FermiDirac(E,chemical_potential,T)* DOS_Gauss_item.N .* GaussDOS(E,DOS_Gauss_item.E_ref,DOS_Gauss_item.param); 

end

