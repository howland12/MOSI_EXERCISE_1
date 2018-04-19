function [ DOS_admin ] = AddGaussToDOS( DOS_admin,E, N_G, E1, E2, ch_type )
% provides Gaussian-shaped DOS  
% in units of states per energy intervall eV per unit volume in m
% provides vector to store the DOS

% input: DOS_admin .. DOS administration data
%        N_G .. density of states
%        E   .. energy vector
%        E1  .. energy of Gaussian center = mean value (in eV)
%        E2  .. width of Gaussian = "sigma"            (in eV)
%        type .. Donor or acceptor?

% output: DOS vector with number of states per energy interval 
%         corresponding to requested Gaussian added
%         occupation vector (zero-valued)

% check whether donor or acceptor like
switch ch_type
    case 'P'    % acceptor
        label = 'AG';
    case 'N'    % donor
        label = 'DG';
    case '0'    % neutral
        label = '0G';
end 

    DOS_temp = zeros(length(E),'like',E);
    DOS_temp = N_G * 1/sqrt(2*pi)/E2* exp(-((E-E1)/E2).^2/2 );

DOS_admin = AddContribtionToDOS(DOS_admin, DOS_temp, label,ch_type,E1,E2,N_G);

%occ_vector = zeros(size(E),'like',E);
end



