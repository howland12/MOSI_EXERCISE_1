%
%   this function emulates the charge neutrality condition
%   desired shape F(E) with F(E) = 0 when condition holds
%   -> one evaluates the expression F(E) = n(E)-p(E)

% input:    Energy    .. reference energy at which condition ought to be
%                        evaluated  / eV
%           DOS_admin .. complete information on DOS
%           m_n_eff   .. effective electron mass
%           m_p_eff   .. effective hole mass
%           T         .. temperature / K
%           

function [ eval ] = chargeNeutrality(Energy, DOS_admin, m_n_eff,m_p_eff,T)

    N = 0;      % total density associated to negative charges / inverse cubic meters
    P = 0;      % total density associated to posiitve charges / inverse cubic meters

% here we need to collect all contributions relevant for charge neutrality

    number_of_DOS_entries = size(DOS_admin,2);

% check for type of DOS contribution

    for k=1:number_of_DOS_entries
        
        n = 0; p = 0; % temporary densities
        
        switch (DOS_admin(k).label(2))
    
            case 'B'  % band state (sqrt-shaped)

                if (DOS_admin(k).label(1) == 'C')
                    n = GetDensityInBand(Energy,DOS_admin(k).E_ref, m_n_eff,T);
                elseif (DOS_admin(k).label(1) == 'V')
                    p = GetDensityInBand(Energy,DOS_admin(k).E_ref, m_p_eff,T);
                end;

            case 'L'    % sharp level (delta-shaped)

                if (DOS_admin(k).type == 'N')     % acceptor-like
                    n = GetDensityInLevel(Energy,DOS_admin(k),T);
                elseif (DOS_admin(k).type == 'P') % donor-like
                    p = DOS_admin(k).N - GetDensityInLevel(Energy,DOS_admin(k),T);
                end
              
            case 'G'    % Gaussian-shaped

                if (DOS_admin(k).type == 'N')   % acceptor-like
                    n = GetDensityInGauss(Energy,DOS_admin(k),T);
                elseif (DOS_admin(k).type == 'P') % donor-like
                    p = DOS_admin(k).N - GetDensityInGauss(Energy,DOS_admin(k),T);
                end

        end; % switch

        N = N + n; 
        P = P + p;
    end;

% rather than checking for the absolute deviation from zero
% look into relative error 

    if (N+P ~= 0)
        eval = (N-P)/(N+P);
    else
        eval = 0.0;
    end;
end

