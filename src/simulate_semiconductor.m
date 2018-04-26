function [struct_semiconductor] = simulate_semiconductor( semiconductor, doping_energy, dopant_density, doping_type)


    % -----------------------------------------------------------------------
    %  determination of Fermi levels and charge carrier densities
    % -----------------------------------------------------------------------
    %
    % hints for usage
    % 
    % + all spatial quantities are given in SI
    % + temperatures are given in K
    % + all energies are given in eV

    %-------------------------------------------------------------------------
    %--- parameter section ---------------------------------------------------
    %-------------------------------------------------------------------------
    % 
    k_eV = 1.3806504E-23 / 1.602176487E-19; % Boltzmann constant in eV/K

    % external conditions
    T = 300;  % in K

    % technical parameters
    % number of intervals in DOS and energy
    resolution = 1500;

    % threshold and maximal number of iterations for root finding
    tolerance = 1d-4;
    max_RF_iter = 30;

    % script control
    plotlog = true;
    SetPlotProperties() ;
    %-------------------------------------------------------------------------
    %--- begin script --------------------------------------------------------
    %-------------------------------------------------------------------------


    % (1) Initialize vectors storing variations

    % temperature 0 .. 800 K
    temperature = linspace(10,800,100); % vector with temperatures in K

    % initalize vector storing electron densities
    main_charge_carrier_number = zeros(size(temperature),'like',temperature);
    n_i = zeros(size(temperature),'like',temperature);
    ionized_dopants = zeros(size(temperature),'like',temperature);

    chemical_potential = zeros(size(temperature),'like',temperature);
    chemical_potential_i = zeros(size(temperature),'like',temperature);

    %-------------------------------------------------------------------------
    % (2a) Initialize energy interval, DOS, and occupation vector 
    %      -> available as column vectors
    %-------------------------------------------------------------------------


    [E_C, E_V, m_n_eff, m_p_eff] = AssignSemiconductor(semiconductor);

    % size of energy interval / eV 
    E_min = E_V - 0.5;
    E_max = E_C + 0.5 ;

    [energies,DOS, occupation] = InitializeEnergyAndDOS(E_min, E_max, E_V, ...
                                                        E_C,resolution);

    DOS_admin = InitializeDOSAdministration(energies);

    % add conduction band, add empty vector for electron occupation
    DOS_admin = AddConductionBandToDOS(DOS_admin,energies,E_C,m_n_eff);

    % add hole band, add empty vector for hole occupation
    DOS_admin = AddValenceBandToDOS(DOS_admin,energies,E_V,m_p_eff);

    % add donor level 100meV above valence band, 
    % level become positive upon emptying
    
    if doping_type == 'p-doped'
        DOS_admin  = AddLevelToDOS(DOS_admin,energies,dopant_density,...
                               E_V + doping_energy, 'N');
    else
        DOS_admin  = AddLevelToDOS(DOS_admin,energies,dopant_density,...
                               E_C + doping_energy, 'P');
    end


    %-------------------------------------------------------------------------
    % (4a) investigate impact of temperature
    %-------------------------------------------------------------------------


    num_temperatures = length(temperature);
    % for each temperature

    for k=1:num_temperatures

        % --------------------------------------------------------------------
        %  evaluate intrinsic Fermi level numerically
        % --------------------------------------------------------------------
        % (a) cast charge neutrality condition into a form F(E,...) = 0 
        % (b) function F has to be provided, here chargeNeutralityIntrinsic()
        % (c) pass F as function handle fh, make sure that E is indicated as the
        %     argument to be evaluated

        fh = @(E) chargeNeutralityIntrinsic(E ,E_C,E_V,m_n_eff,m_p_eff,...
                                            temperature(k));

        % (d) employ root-finding algorithm to determine the chemical potential

        [chemical_potential_i(k), num_iter, error] = ...
            FindRootNestedIntervals(fh,energies, (E_C + E_V)/2.+0.2,...
            tolerance, max_RF_iter);

        n = GetDensityInBand(chemical_potential_i(k),E_C,m_n_eff, ...
                             temperature(k));
        p = GetDensityInBand(chemical_potential_i(k),E_V,m_p_eff, ...
                             temperature(k));
        n_i(k) = sqrt(n*p);

        %---------------------------------------------------------------------
        % evaluate Fermi level numerically for a non-intrinsic system
        %---------------------------------------------------------------------

        % (a) cast charge neutrality condition into a form F(E,...) = 0 
        % (b) function F has to be provided, here chargeNeutrality()
        % (c) pass F as function handle fh, make sure that E is indicated as the
        %     argument to be evaluated

        fh = @(E) chargeNeutrality(E,DOS_admin,m_n_eff,m_p_eff,temperature(k));

        % (d) employ root-finding algorithm to determine the chemical potential

        [chemical_potential(k), num_iter, error] = ...
            FindRootNestedIntervals(fh,energies, ...
            chemical_potential_i(k), tolerance, max_RF_iter);

        main_charge_carrier_number(k) = GetDensityInBand(chemical_potential(k), ...
                                     E_C,m_n_eff, temperature(k));
        if doping_type == 'p-doped'
            
            main_charge_carrier_number(k) = GetDensityInBand(chemical_potential(k), ...
                                     E_V, m_p_eff, temperature(k));
            
            ionized_dopants(k) = GetDensityInLevel(chemical_potential(k),...
                                 DOS_admin(3),temperature(k))/DOS_admin(3).N;
                             
        else
            
            main_charge_carrier_number(k) = GetDensityInBand(chemical_potential(k), ...
                                     E_C,m_n_eff, temperature(k));
            
            ionized_dopants(k) = 1.0 - ...
                                 GetDensityInLevel(chemical_potential(k),...
                                 DOS_admin(3),temperature(k))/DOS_admin(3).N;
                             
        end
        

    end;
         
    
    struct_semiconductor.E_V = E_V;
    struct_semiconductor.E_C = E_C;
    struct_semiconductor.temperature = temperature;
    struct_semiconductor.dopant_density = dopant_density;
    struct_semiconductor.n_i = n_i;
    struct_semiconductor.chemical_potential_i = chemical_potential_i;
    struct_semiconductor.chemical_potential = chemical_potential;
    struct_semiconductor.ionized_dopants = ionized_dopants;
    struct_semiconductor.main_charge_carrier_number = main_charge_carrier_number;
    

    %%

