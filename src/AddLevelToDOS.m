function [ DOS_admin ] = AddLevelToDOS( DOS_admin,E, N_G, E_level, ch_type )
% provides Gaussian-shaped DOS  
% in units of states per energy intervall eV per unit volume in m
% provides vector to store the DOS

% input: DOS_admin .. DOS administration data
%        N_G .. density of states
%        E   .. energy vector
%        E_level  .. energy of Gaussian center = mean value (in eV)
%        type .. Donor or acceptor?

% output: DOS vector with number of states per energy interval 
%         corresponding to requested Gaussian added
%         occupation vector (zero-valued)

% check whether donor or acceptor like
switch ch_type
    case 'P'    % donor
        label = 'DL';
    case 'N'    % acceptor
        label = 'AL';
    case '0'    % neutral
        label = '0L';
end 

    DOS_temp = zeros(size(E),'like',E);

    % add delta-shaped level
    interval_index = max(find(E<E_level));
    DOS_temp(interval_index) = N_G;

    DOS_admin = AddContribtionToDOS(DOS_admin, DOS_temp, label,ch_type,E_level,0,N_G);

%occ_vector = zeros(size(E),'like',E);
end



