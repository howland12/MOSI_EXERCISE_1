
clearvars;
close all;

% temperature 0 .. 800 K
temperature = linspace(10,800,100); % vector with temperatures in K

[Si_p_doped] = simulate_semiconductor('Si', 0.1, 'sharp', 1e21,'p-doped', temperature);
[Ge_p_doped] = simulate_semiconductor('Ge', 0.1, 'sharp', 1e21,'p-doped', temperature);
[GaAs_p_doped] = simulate_semiconductor('GaAs', 0.1, 'sharp', 1e21,'p-doped', temperature);

clf(figure(1))
h=figure(1);
         
         set(gcf,'color','w');
        hold on

        plot(Si_p_doped.temperature,Si_p_doped.E_V * ones(size(Si_p_doped.temperature)),... 
             'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_V');
        plot(Si_p_doped.temperature,Si_p_doped.E_C * ones(size(Si_p_doped.temperature)),... 
             'LineWidth',1,'Color','k','DisplayName','Si E_C');
        plot(Si_p_doped.temperature,Si_p_doped.chemical_potential_i,'--','LineWidth',1,...
             'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
        plot(Si_p_doped.temperature,...
            Si_p_doped.chemical_potential,...
            'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');

        title({'chemical potential vs temperature',' ',...
               'in Si at N_A = 10^{21} m^3'});
           ax = gca;
ax.FontSize = 11;
grid on
        legend('E_V','E_C','\mu_i','\mu', 'Location' ,'northeastoutside');
 
        ylim([0 1.4]);
        xlabel('temperature / K');
        ylabel('energy / eV');
        
% figure 
clf(figure(2))
figure(2)

    hold on

    plot(Si_p_doped.temperature,  Si_p_doped.n_i/Si_p_doped.dopant_density,'.','LineWidth',1,'Color',...
         [1 0 0],'DisplayName','Si_p_doped');
    plot(Si_p_doped.temperature, Si_p_doped.main_charge_carrier_number/Si_p_doped.dopant_density,'LineWidth',2,'Color',...
         [1 0 0],'DisplayName','Si_p_doped');

    title({'electron density vs temperature',' ',...
          'in GaAs at N_D = 10^{21} m^3'});
                ax = gca;
ax.FontSize = 11;
grid on
    legend('n_i/N_D','n/N_D', 'Location' ,'northeastoutside');

    ylim([0 2.0]);
    xlabel('temperature / K');
    ylabel('electron density / N_D');

        
        
        
        
% clf(figure(2))
% figure(2)
%         hold on
% 
%         plot(Ge_p_doped.temperature,Ge_p_doped.E_V * ones(size(Ge_p_doped.temperature)),... 
%              'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_C');
%         plot(Ge_p_doped.temperature,Ge_p_doped.chemical_potential_i,'--','LineWidth',1,...
%              'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
%         plot(Ge_p_doped.temperature(find(Ge_p_doped.chemical_potential > Ge_p_doped.E_V)),...
%             Ge_p_doped.chemical_potential(find(Ge_p_doped.chemical_potential > Ge_p_doped.E_V)),...
%             'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');
% 
%         title({'chemical potential vs temperature',' ',...
%                'in Si at N_A = 10^{21} m^3'});
%         legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');
% 
%         ylim([0 1.5]);
%         xlabel('temperature / K');
%         ylabel('energy / eV');
%         
% clf(figure(3))
% figure(3)
%         hold on
% 
%         plot(GaAs_p_doped.temperature,GaAs_p_doped.E_V * ones(size(GaAs_p_doped.temperature)),... 
%              'LineWidth',1,'Color',[0 0 1],'DisplayName','Si E_C');
%         plot(GaAs_p_doped.temperature,GaAs_p_doped.chemical_potential_i,'--','LineWidth',1,...
%              'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
%         plot(GaAs_p_doped.temperature(find(GaAs_p_doped.chemical_potential > GaAs_p_doped.E_V)),...
%             GaAs_p_doped.chemical_potential(find(GaAs_p_doped.chemical_potential > GaAs_p_doped.E_V)),...
%             'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');
% 
%         title({'chemical potential vs temperature',' ',...
%                'in Si at N_A = 10^{21} m^3'});
%         legend('E_C','\mu_i','\mu', 'Location' ,'northeastoutside');
% 
%         ylim([0 1.5]);
%         xlabel('temperature / K');
%         ylabel('energy / eV');
% 
% 
% clf(figure(4))    
% figure(4)
% 
%     hold on
% 
%     plot(Si_p_doped.temperature, Si_p_doped.n_i/Si_p_doped.dopant_density,'.','LineWidth',1,'Color',...
%          [1 0 0],'DisplayName','Si');
%      plot(Ge_p_doped.temperature, Ge_p_doped.n_i/Ge_p_doped.dopant_density,'.','LineWidth',1,'Color',...
%          'b','DisplayName','Si');
%      plot(GaAs_p_doped.temperature, GaAs_p_doped.n_i/GaAs_p_doped.dopant_density,'.','LineWidth',1,'Color',...
%          'g','DisplayName','Si');
%     plot(Si_p_doped.temperature, Si_p_doped.main_charge_carrier_number/Si_p_doped.dopant_density,'LineWidth',2,'Color',...
%          [1 0 0],'DisplayName','Si');
% 
%     title({'hole density vs temperature',' ',...
%           'in Si at N_A = 10^{21} m^3'}); 
%     legend('n_i/N_A','p/N_A', 'Location' ,'northeastoutside');  
% 
%     ylim([0 2.0]);
%     xlabel('temperature / K');
%     ylabel('hole density / N_A');
% 
% clf(figure(5))    
% figure(5)
% 
% 
%     % correct for evalation errors
%     Si_p_doped.ionized_dopants(find(Si_p_doped.ionized_dopants == 1)) = 0;
% 
%     hold on
%     title({'number of ionized dopants vs temperature',' ',...
%            'in Si at N_A = 10^{21} m^3'}); 
% 
%     plot(Si_p_doped.temperature, Si_p_doped.ionized_dopants,'LineWidth',2,'Color',[1 0 0],...
%          'DisplayName','Ge');
%     legend('N_A^+/N_A', 'Location' ,'northeastoutside');
% 
%     ylim([0 1.1]);
%     xlabel('temperature / K');
%     ylabel('density of ionized dopants/ N_A');
%     
%   
%     
%     
% [Si_n_doped] = simulate_semiconductor('Si', -0.1, 1e21,'n-doped');
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         TASK 2                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- common parameters for task 2 -----------------------
temperature = 300;
ND_vector = logspace(13,25,100);


%############################## a #########################################

%-------------------------- simulation loop -> task a) --------------------
data_container_a = cell(1,length(ND_vector));
for i = 1 : length(ND_vector)
    data_container_a{1,i} = simulate_semiconductor('Si', -0.660000000000000, 'sharp', ND_vector(i),'n-doped', temperature);
end

%---------------------- extract needed data from container -> task a) -----
chemical_potential_vector_a = zeros([1 length(ND_vector)]);
for i = 1 : length(ND_vector)
    chemical_potential_vector_a(i) = data_container_a{i}.chemical_potential;
end

%############################## b #########################################

%-------------------------- simulation loop -> task b) --------------------
donor_state_energy = -(data_container_a{1}.E_C / 2 + 0.1);
ND_vector = logspace(13,25,100);
data_container_b = cell(1,length(ND_vector));
for i = 1 : length(ND_vector)
    data_container_b{1,i} = simulate_semiconductor('Si', donor_state_energy, 'gaussian', ND_vector(i),'n-doped', temperature);
end

%---------------------- extract needed data from container -> task b) -----
chemical_potential_vector_b = zeros([1 length(ND_vector)]);
for i = 1 : length(ND_vector)
    chemical_potential_vector_b(i) = data_container_b{i}.chemical_potential;
end

%###################### Results of Task 2 a) and b) ####################### 

%------------------------- plot simulation results -> task a) --------------
figure();
set(gcf,'color','w');
figure()
semilogx(ND_vector, chemical_potential_vector_a);
grid on
xlabel('N_D / 1');
ylabel('\mu / ev');
xlim([ND_vector(1), ND_vector(end)]);

%------------------------- plot simulation results -> task b) -------------
figure();
set(gcf,'color','w');
figure()
semilogx(ND_vector, chemical_potential_vector_b);
grid on
xlabel('N_D / 1');
ylabel('\mu / ev');
xlim([ND_vector(1), ND_vector(end)]);


%------------------------ plot combined results -> task a) and b) ---------
figure();
semilogx(ND_vector, chemical_potential_vector_a);
hold on;
semilogx(ND_vector, chemical_potential_vector_b);
hold off;
grid on
xlabel('N_D / 1');
ylabel('\mu / ev');
xlim([ND_vector(1), ND_vector(end)]);
legend('sharp level shape','gaussian shape')
    
    
    
    
    
    
    
    