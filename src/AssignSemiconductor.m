function [E_C, E_V, m_n_eff, m_p_eff] = AssignSemiconductor(label);

% intrinsic semiconductor properties


% input:   label of compound (4 character string)
% output:  E_C     .. conduction band minimum / eV
%          E_V     .. valence band minimum / eV
%          m_n_eff .. effective electron mass / m_e
%          m_h_eff .. effective hole mass / m_e

%
% for a collection of material parameters, see
% http://lampx.tugraz.at/~hadley/psd/L3/L3.php
% https://ecee.colorado.edu/~bart/book/book/append/append3.htm


% within this example we work with the following approximation:
%
% For the valence and conduction band, ONE isotropic parabolic band with
% ONE associated effective mass is assumed
%
% i.e., we disgard
% + that the effective mass depends on the direction of motion, as bands do
%   not necessarily adopt the same curvature along all possible k
%   directions; longitudinal and transversal effective masses are not
%   discriminated
%   ml* = 1.64  ;  mt* = 0.082
%
% + the presence of multiple bands that share the same maximum or minimum
%   energy, e.g., the heavy and light hole bands of silicon
%   Si: mlh* = 0.044;  mhh* = 0.28 
%
% Consequently, the values of effective masses given below describe best a
% square-root-shaped DOS containing all contributing bands and directions


% list of effective masses (in units of electron mass)
%   element          Ge           Si           GaAs
%   -------------+--------+-------------+-------------
%   electron        0.55        1.08           0.067  (Si = 1.18?)
%   hole            0.37        0.81           0.45
%   hole (heavy)


% list of band gaps (in eV)
%   element          Ge           Si           GaAs
%   -------------+--------+-------------+-------------
%   band gap        0.66         1.12           1.424


E_V = 0. ;  % eV
E_C = 0. ;
m_n_eff = 1.;
m_p_eff = 1.;

    while (length(label) < 4)
        label = [' ',label];
    end;

    switch label

        case '  Ge'  % Germanium

            m_n_eff = 0.55 ;    % in units of electron mass
            m_p_eff = 0.37 ;    % in units of electron mass
            E_C = E_V + 0.66;   % eV

        case '  Si'  % Si

            m_n_eff = 1.18 ;    % in units of electron mass
            m_p_eff = 0.81 ;    % in units of electron mass
            E_C = E_V + 1.12;   % eV

        case 'GaAs'  % GaAs

            m_n_eff = 0.067 ;   % in units of electron mass
            m_p_eff = 0.45 ;    % in units of electron mass
            E_C = E_V + 1.424;  % eV
            
        case 'In3As'  % In3As: http://www.ioffe.ru/SVA/NSM/Semicond/InAs/bandstr.html
            
            m_n_eff = 0.023 ;   % in units of electron mass
            m_p_eff = 0.41 ;    % in units of electron mass
            E_C = E_V + 0.354;  % eV

    end; %switch
end