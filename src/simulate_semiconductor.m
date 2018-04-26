function [E_V, E_C, temperature, dopant_density, n_i, n, p, chemical_potential_i, chemical_potential, ionized_dopants, main_charge_carrier_number] = simulate_semiconductor( semiconductor, doping_energy, doping_type)


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

    dopant_density = 1e21; % in SI

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

    %%


 

% clf(figure(7))
% figure(7)
%         hold on
% 
%         plot(temperature,E_V * ones(size(temperature)),... 
%              'LineWidth',1,'Color',[0 0 1],'DisplayName','Ge E_C');
%         plot(temperature,chemical_potential_i,'--','LineWidth',1,...
%              'Color',[1 0 0],'DisplayName','Ge E_F_intrinsic');
%         plot(temperature(find(chemical_potential < E_C)),...
%             chemical_potential(find(chemical_potential < E_C)),...
%             'LineWidth',2,'Color',[0 1 0],'DisplayName','Ge \mu');
% 
%         title({'chemical potential vs temperature',' ',...
%                'in Ge at N_D = 10^{21} m^3'});
%         legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');
% 
%         ylim([0 1.5]);
%         xlabel('temperature / K');
%         ylabel('energy / eV');


%     clf(figure(8))    
%     figure(8)
% 
%         hold on
% 
%         plot(temperature, n_i/dopant_density,'.','LineWidth',1,'Color',...
%              [1 0 0],'DisplayName','Ge');
%         plot(temperature, main_charge_carrier_number/dopant_density,'LineWidth',2,'Color',...
%              [1 0 0],'DisplayName','Ge');
% 
%         title({'electron density vs temperature',' ',...
%               'in Ge at N_D = 10^{21} m^3'}); 
%         legend('n_i/N_D','n/N_D', 'Location' ,'northeastoutside');  
% 
%         ylim([0 2.0]);
%         xlabel('temperature / K');
%         ylabel('electron density / N_D');


    clf(figure(9))    
    figure(9)


        % correct for evalation errors
        ionized_dopants(find(ionized_dopants == 1)) = 0;

        hold on
        title({'number of ionized dopants vs temperature',' ',...
               'in Ge at N_D = 10^{21} m^3'}); 

        plot(temperature, ionized_dopants,'LineWidth',2,'Color',[1 0 0],...
             'DisplayName','Ge');
        legend('N_D^+/N_D', 'Location' ,'northeastoutside');

        ylim([0 1.1]);
        xlabel('temperature / K');
        ylabel('density of ionized dopants/ N_D');

