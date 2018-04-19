function [ Fhalf ] = FintegrationForBands( eta )
    %func_handle = @(E)integrandForBands(E, E_offset);

% upper boundary for integration is infinity
% approximate this boundary with a large value
    % y = integration(func_handle,0,1E4)

% use here approximate analytical expression taken from 
% https://www.eecis.udel.edu/~kolodzey/courses/ELEG667F06/667F06HW/Pierret_6_FermiDirac.pdf
% Table 4.2 (p. 118)

a = (power(abs(eta-2.13),2.4) + 9.6);
a = power(a,5/12) + eta + 2.13;
a = power(a,-1.5)* 3*sqrt(pi/2) + exp(-eta);


% suggestion of WIKIPEDIA.org
%
%if (eta < 1.3)
%    a = exp(-eta) + 2.7;
%end;    
   
%


Fhalf = 1/a;
end

