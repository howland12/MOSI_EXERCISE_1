function [ eval ] = chargeNeutralityIntrinsic(Energy, E_C, E_V, m_n_eff,m_p_eff,T)
%   this function emulates the charge neutrality condition
%   desired shape F(E) with F(E) = 0 when condition holds
%   -> one evaluates the expression F(E) = n(E)-p(E)

    n = GetDensityInBand(Energy,E_C, m_n_eff,T);

    p = GetDensityInBand(Energy,E_V, m_p_eff,T);

% rather than checking for the absolute deviation from zero
% look into relative error 

    eval = (n-p)/(n+p);
end

