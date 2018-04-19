% gives DOS at a specific energy for a GAUSSIAN-shaped DOS

% input: E .. energy
%        Ecenter .. shape-defining parameters, here mean value
%        Ewidth  .. shape-defining parameter, here width sigma


function [ dos_val ] = GaussDOS( E, Ecenter, Ewidth )


dos_val = 1/sqrt(2*pi)/Ewidth* exp(  -((E-Ecenter)/Ewidth).^2/2 );


end

