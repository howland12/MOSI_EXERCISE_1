function [ prob ] = FermiDirac( E, CP,T )
% Gives occupation probability with FERMI DIRAC distribution
% input: E   .. energy
%        CP .. chemical potential
%        T   .. temperature
% output: prob .. occupation probability

k = 1.3806504d-23; %in SI
q = 1.602176565e-19; 

% convert eV to SI vy using a temperature equivalent
T = T/q;

prob = 1./(exp((E-CP)/k/T)+1);


end

