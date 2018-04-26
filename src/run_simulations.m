close all;
clearvars;

run_GaAs = true;
run_Si = true;
run_Ge = true;

% if run_GaAs
%     run('main_GaAs.m')
% end
% 
% if run_Si
%     run('main_Si.m')
% end
% 
% if run_Ge
%     run('main_Ge.m')
% end

[E_V, E_C, temperature, dopant_density, n_i, n, p, chemical_potential_i, chemical_potential, ionized_dopants, main_charge_carrier_number] = simulate_semiconductor('Si', 0.1, 'p-doped');
clf(figure(7))
figure(7)
        hold on

        plot(temperature,E_V * ones(size(temperature)),... 
             'LineWidth',1,'Color',[0 0 1],'DisplayName','Ge E_C');
        plot(temperature,chemical_potential_i,'--','LineWidth',1,...
             'Color',[1 0 0],'DisplayName','Ge E_F_intrinsic');
        plot(temperature(find(chemical_potential < E_C)),...
            chemical_potential(find(chemical_potential < E_C)),...
            'LineWidth',2,'Color',[0 1 0],'DisplayName','Ge \mu');

        title({'chemical potential vs temperature',' ',...
               'in Ge at N_D = 10^{21} m^3'});
        legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

        ylim([0 1.5]);
        xlabel('temperature / K');
        ylabel('energy / eV');


clf(figure(8))    
figure(8)

    hold on

    plot(temperature, n_i/dopant_density,'.','LineWidth',1,'Color',...
         [1 0 0],'DisplayName','Ge');
    plot(temperature, main_charge_carrier_number/dopant_density,'LineWidth',2,'Color',...
         [1 0 0],'DisplayName','Ge');

    title({'electron density vs temperature',' ',...
          'in Ge at N_D = 10^{21} m^3'}); 
    legend('n_i/N_D','n/N_D', 'Location' ,'northeastoutside');  

    ylim([0 2.0]);
    xlabel('temperature / K');
    ylabel('electron density / N_D');

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