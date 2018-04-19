
% -----------------------------------------------------------------------
%  determination of Fermi levels and charge carrier densities
% -----------------------------------------------------------------------
%
% hints for usage
% 
% + all spatial quantities are given in SI
% + temperatures are given in K
% + all energies are given in eV

% clear figures
% close all
% clear all
clearvars -except run_GaAs run_Si run_Ge

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
n_GaAs = zeros(size(temperature),'like',temperature);
n_i_GaAs = zeros(size(temperature),'like',temperature);
ND_ionized_GaAs = zeros(size(temperature),'like',temperature);

chemical_potential_GaAs = zeros(size(temperature),'like',temperature);
chemical_potential_i_GaAs = zeros(size(temperature),'like',temperature);

%-------------------------------------------------------------------------
% (2a) Initialize energy interval, DOS, and occupation vector 
%      -> available as column vectors
%-------------------------------------------------------------------------


[E_C, E_V, m_n_eff, m_p_eff] = AssignSemiconductor('GaAs');

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

% add donor level 50meV below conduction band, 
% level become positive upon emptying
DOS_admin  = AddLevelToDOS(DOS_admin,energies,dopant_density,...
                           E_C - 0.05, 'P');

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

    [chemical_potential_i_GaAs(k), num_iter, error] = ...
        FindRootNestedIntervals(fh,energies, (E_C + E_V)/2.+0.2,...
        tolerance, max_RF_iter);

    n = GetDensityInBand(chemical_potential_i_GaAs(k),E_C,m_n_eff, ...
                         temperature(k));
    p = GetDensityInBand(chemical_potential_i_GaAs(k),E_V,m_p_eff, ...
                         temperature(k));
    n_i_GaAs(k) = sqrt(n*p);

    %---------------------------------------------------------------------
    % evaluate Fermi level numerically for a non-intrinsic system
    %---------------------------------------------------------------------

    % (a) cast charge neutrality condition into a form F(E,...) = 0 
    % (b) function F has to be provided, here chargeNeutrality()
    % (c) pass F as function handle fh, make sure that E is indicated as the
    %     argument to be evaluated

    fh = @(E) chargeNeutrality(E,DOS_admin,m_n_eff,m_p_eff,temperature(k));

    % (d) employ root-finding algorithm to determine the chemical potential

    [chemical_potential_GaAs(k), num_iter, error] = ...
        FindRootNestedIntervals(fh,energies, ...
        chemical_potential_i_GaAs(k), tolerance, max_RF_iter);

    n_GaAs(k) = GetDensityInBand(chemical_potential_GaAs(k), ...
                                 E_C,m_n_eff, temperature(k));

    ND_ionized_GaAs(k) = 1.0 - ...
                         GetDensityInLevel(chemical_potential_GaAs(k),...
                         DOS_admin(3),temperature(k))/DOS_admin(3).N;

end;

%%


clf(figure(1))
figure(1)
    
    plot(temperature,E_C * ones(size(temperature),'like',dopant_density),... 
         'LineWidth',1,'Color',[0 1 0],'DisplayName','GaAs E_C');
    

    hold on;
    
    plot(temperature,chemical_potential_i_GaAs,'--','LineWidth',1,...
         'Color',[0 0 1],'DisplayName','GaAs E_F_intrinsic');
     
    plot(temperature(find(chemical_potential_GaAs < E_C)),...
        chemical_potential_GaAs(find(chemical_potential_GaAs < E_C)),...
        'LineWidth',2,'Color',[1 0 0],'DisplayName','GaAs \mu');
    
    hold off;
    
    
    
    title({'chemical potential vs temperature',' ',...
           'in GaAs at N_D = 10^{21} m^3'});
       
    legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');
    
    ylim([0 1.5]);
    xlabel('temperature / K');
    ylabel('energy / eV');
     

    
clf(figure(2))    
figure(2)

    hold on
    
    plot(temperature, n_i_GaAs/dopant_density,'.','LineWidth',1,'Color',...
         [1 0 0],'DisplayName','GaAs');
    plot(temperature, n_GaAs/dopant_density,'LineWidth',2,'Color',...
         [1 0 0],'DisplayName','GaAs');
    
    title({'electron density vs temperature',' ',...
          'in GaAs at N_D = 10^{21} m^3'}); 
    legend('n_i/N_D','n/N_D', 'Location' ,'northeastoutside');  
           
    ylim([0 2.0]);
    xlabel('temperature / K');
    ylabel('electron density / N_D');

    
clf(figure(3))    
figure(3)


    % correct for evalation errors
    ND_ionized_GaAs(find(ND_ionized_GaAs == 1)) = 0;
    
    hold on
    title({'number of ionized dopants vs temperature',' ',...
           'in GaAs at N_D = 10^{21} m^3'}); 

    plot(temperature, ND_ionized_GaAs,'LineWidth',2,'Color',[1 0 0],...
         'DisplayName','GaAs');
    legend('N_D^+/N_D', 'Location' ,'northeastoutside');

    ylim([0 1.1]);
    xlabel('temperature / K');
    ylabel('density of ionized dopants/ N_D');
    
