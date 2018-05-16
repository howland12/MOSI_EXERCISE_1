
clearvars;
close all;

% temperature 0 .. 800 K
temperature = linspace(10,800,100); % vector with temperatures in K

[Si_p_doped] = simulate_semiconductor('Si', 0.1, 1e21,'p-doped', temperature);
[Ge_p_doped] = simulate_semiconductor('Ge', 0.1, 1e21,'p-doped', temperature);
[GaAs_p_doped] = simulate_semiconductor('GaAs', 0.1, 1e21,'p-doped', temperature);

clf(figure(1))
figure(1)
        hold on

        plot(Si_p_doped.temperature,Si_p_doped.E_V * ones(size(Si_p_doped.temperature)),... 
             'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_C');
        plot(Si_p_doped.temperature,Si_p_doped.chemical_potential_i,'--','LineWidth',1,...
             'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
        plot(Si_p_doped.temperature(find(Si_p_doped.chemical_potential > Si_p_doped.E_V)),...
            Si_p_doped.chemical_potential(find(Si_p_doped.chemical_potential > Si_p_doped.E_V)),...
            'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');

        title({'chemical potential vs temperature',' ',...
               'in Si at N_A = 10^{21} m^3'});
        legend('E_V','\mu_i','\mu', 'Location' ,'northeastoutside');

        ylim([0 1.5]);
        xlabel('temperature / K');
        ylabel('energy / eV');
        
clf(figure(2))
figure(2)
        hold on

        plot(Ge_p_doped.temperature,Ge_p_doped.E_V * ones(size(Ge_p_doped.temperature)),... 
             'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_C');
        plot(Ge_p_doped.temperature,Ge_p_doped.chemical_potential_i,'--','LineWidth',1,...
             'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
        plot(Ge_p_doped.temperature(find(Ge_p_doped.chemical_potential > Ge_p_doped.E_V)),...
            Ge_p_doped.chemical_potential(find(Ge_p_doped.chemical_potential > Ge_p_doped.E_V)),...
            'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');

        title({'chemical potential vs temperature',' ',...
               'in Si at N_A = 10^{21} m^3'});
        legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

        ylim([0 1.5]);
        xlabel('temperature / K');
        ylabel('energy / eV');
        
clf(figure(3))
figure(3)
        hold on

        plot(GaAs_p_doped.temperature,GaAs_p_doped.E_V * ones(size(GaAs_p_doped.temperature)),... 
             'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_C');
        plot(GaAs_p_doped.temperature,GaAs_p_doped.chemical_potential_i,'--','LineWidth',1,...
             'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
        plot(GaAs_p_doped.temperature(find(GaAs_p_doped.chemical_potential > GaAs_p_doped.E_V)),...
            GaAs_p_doped.chemical_potential(find(GaAs_p_doped.chemical_potential > GaAs_p_doped.E_V)),...
            'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');

        title({'chemical potential vs temperature',' ',...
               'in Si at N_A = 10^{21} m^3'});
        legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

        ylim([0 0.3]);
        xlabel('temperature / K');
        ylabel('energy / eV');


clf(figure(4))    
figure(4)

    hold on

    plot(Si_p_doped.temperature, Si_p_doped.n_i/Si_p_doped.dopant_density,'.','LineWidth',1,'Color',...
         [1 0 0],'DisplayName','Si');
    plot(Si_p_doped.temperature, Si_p_doped.main_charge_carrier_number/Si_p_doped.dopant_density,'LineWidth',2,'Color',...
         [1 0 0],'DisplayName','Si');

    title({'hole density vs temperature',' ',...
          'in Si at N_A = 10^{21} m^3'}); 
    legend('n_i/N_A','p/N_A', 'Location' ,'northeastoutside');  

    ylim([0 2.0]);
    xlabel('temperature / K');
    ylabel('hole density / N_A');

clf(figure(5))    
figure(5)


    % correct for evalation errors
    Si_p_doped.ionized_dopants(find(Si_p_doped.ionized_dopants == 1)) = 0;

    hold on
    title({'number of ionized dopants vs temperature',' ',...
           'in Si at N_A = 10^{21} m^3'}); 

    plot(Si_p_doped.temperature, Si_p_doped.ionized_dopants,'LineWidth',2,'Color',[1 0 0],...
         'DisplayName','Ge');
    legend('N_A^+/N_A', 'Location' ,'northeastoutside');

    ylim([0 1.1]);
    xlabel('temperature / K');
    ylabel('density of ionized dopants/ N_A');
    
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         TASK 2                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

temperature = 300;

%-------------------------- simulation loop -------------------------------
ND_vector = logspace(13,25,100);
data_container = cell(1,length(ND_vector));
for i = 1 : length(ND_vector)
    data_container{1,i} = simulate_semiconductor('Si', -0.1, ND_vector(i),'n-doped', temperature);
end

%-------------------------- extract needed data from container ------------
chemical_potential_vector = zeros([1 length(ND_vector)]);
for i = 1 : length(ND_vector)
    chemical_potential_vector(i) = data_container{i}.chemical_potential;
end

%------------------------- plot simulation results ------------------------
figure()
set(gcf,'color','w');
semilogx(ND_vector, chemical_potential_vector);
grid on
xlabel('N_D / 1');
ylabel('\mu / ev');
xlim([ND_vector(1), ND_vector(end)]);






 
    
    
    
    
    
    
    
    
    
    