
clearvars;
close all;
addpath('matlab2tikz/src');

% temperature 0 .. 800 K
temperature = linspace(10,800,100); % vector with temperatures in K
bandgap_multiplier = 1; % bandgab of selected semiconductor is multiplied with this constant 
e_m_n_multiplier = 1; % effective mass of electrons of selected semiconductor is multiplied with this constant
e_m_p_multiplier = 1;


[Si_p_doped] = simulate_semiconductor('Si', 0.1, 1e21,'p-doped', temperature, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
[Ge_p_doped] = simulate_semiconductor('Ge', 0.1, 1e21,'p-doped', temperature, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
[GaAs_p_doped] = simulate_semiconductor('GaAs', 0.1, 1e21,'p-doped', temperature, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         TASK 1                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%                         TASK 1(a,b)                               %%
clf(figure(1))
h=figure(1);  
set(gcf,'color','w');
hold on

plot(Si_p_doped.temperature,Si_p_doped.E_V * ones(size(Si_p_doped.temperature)),... 
'LineWidth',2,'Color',[0 0 1],'DisplayName','Si E_V');
plot(Si_p_doped.temperature,Si_p_doped.E_C * ones(size(Si_p_doped.temperature)),... 
'LineWidth',2,'Color','k','DisplayName','Si E_C');
plot(Si_p_doped.temperature,Si_p_doped.chemical_potential_i,'--','LineWidth',2,...
'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
plot(Si_p_doped.temperature,...
Si_p_doped.chemical_potential,...
'LineWidth',2,'Color',[0 1 0],'DisplayName','Si \mu');

% title({'chemical potential vs temperature',' ',...
% 'in Si at N_A = 10^{21} m^3'});
ax = gca;
ax.FontSize = 11;
grid on
legend('E_V','E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

ylim([0 1.2]);
xlabel('temperature / K');
ylabel('energy / eV');
matlab2tikz('results/Si_p_doped_N_A_=1e21_chemical_potential_vs_temperature.tex','showInfo', false);


clf(figure(2))
figure(2)

hold on
set(gcf,'color','w');
plot(Si_p_doped.temperature,  Si_p_doped.n_i/Si_p_doped.dopant_density,'.','LineWidth',2,'Color',...
[1 0 0],'DisplayName','Si_p_doped');
plot(Si_p_doped.temperature, Si_p_doped.main_charge_carrier_number/Si_p_doped.dopant_density,'LineWidth',2,'Color',...
[1 0 0],'DisplayName','Si_p_doped');

% title({'',' ',...
% 'in Si-p-doped at N_A = 10^{21} m^3'});
ax = gca;
ax.FontSize = 11;
grid on
legend('n_i/N_A','p/N_A', 'Location' ,'northeastoutside');

ylim([0 2.0]);
xlabel('temperature / K');
ylabel('hole density / N_A');
matlab2tikz('results/Si_p_doped_N_A_=1e21_electron_density_vs_temperature.tex','showInfo', false);

    
clf(figure(3))
figure(3)
set(gcf,'color','w');
% correct for evalation errors
Si_p_doped.ionized_dopants(find(Si_p_doped.ionized_dopants == 1)) = 0;

hold on

plot(Si_p_doped.temperature, Si_p_doped.ionized_dopants,'LineWidth',2,'Color',[1 0 0],...
 'DisplayName','Si');
legend('N_A^-/N_A', 'Location' ,'northeastoutside');

% title({'number of ionized dopants vs temperature',' ',...
%    'in Si-p-doped at N_A = 10^{21} m^3'}); 

ax = gca;
ax.FontSize = 11;
grid on

ylim([0 1.1]);
xlabel('temperature / K');
ylabel('density of ionized dopants/ N_A');
matlab2tikz('results/Si_p_doped_N_A_=1e21_number_of_ionized_dopants_vs_temperature.tex','showInfo', false);
    
%%%                         TASK 1(c)                                      %%%     


% chemical potential doped Ge%       
clf(figure(4))
h=figure(4);

set(gcf,'color','w');
hold on

plot(Ge_p_doped.temperature,Ge_p_doped.E_V * ones(size(Si_p_doped.temperature)),... 
 'LineWidth',2,'Color',[0 0 1],'DisplayName','Ge E_V');
plot(Ge_p_doped.temperature,Ge_p_doped.E_C * ones(size(Si_p_doped.temperature)),... 
 'LineWidth',2,'Color','k','DisplayName','Ge E_C');
plot(Ge_p_doped.temperature,Ge_p_doped.chemical_potential_i,'--','LineWidth',2,...
 'Color',[1 0 0],'DisplayName','Ge E_F_intrinsic');
plot(Ge_p_doped.temperature,...
Ge_p_doped.chemical_potential,...
'LineWidth',2,'Color',[0 1 0],'DisplayName','Ge \mu');

% title({'chemical potential vs temperature',' ',...
%    'in Ge at N_A = 10^{21} m^3'});
ax = gca;
ax.FontSize = 11;
grid on
legend('E_V','E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

ylim([0 0.7]);
xlabel('temperature / K');
ylabel('energy / eV');    
matlab2tikz('results/Ge_p_doped_N_A_=1e21_chemical_potential_vs_temperature.tex','showInfo', false);
    

% chemical potential doped Ge, GaAs%      
clf(figure(5))
h=figure(5);

set(gcf,'color','w');
hold on

plot(GaAs_p_doped.temperature,GaAs_p_doped.E_V * ones(size(Si_p_doped.temperature)),... 
 'LineWidth',2,'Color',[0 0 1],'DisplayName','GaAs E_V');
plot(GaAs_p_doped.temperature,GaAs_p_doped.E_C * ones(size(Si_p_doped.temperature)),... 
 'LineWidth',2,'Color','k','DisplayName','GaAs E_C');
plot(GaAs_p_doped.temperature,GaAs_p_doped.chemical_potential_i,'--','LineWidth',2,...
 'Color',[1 0 0],'DisplayName','GaAs E_F_intrinsic');
plot(GaAs_p_doped.temperature,...
GaAs_p_doped.chemical_potential,...
'LineWidth',2,'Color',[0 1 0],'DisplayName','GaAs \mu');

% title({'chemical potential vs temperature',' ',...
%    'in GaAs at N_A = 10^{21} m^3'});
ax = gca;
ax.FontSize = 11;
grid on
legend('E_V','E_C','\mu_i','\mu', 'Location' ,'northeastoutside');

ylim([0 1.5]);
xlabel('temperature / K');
ylabel('energy / eV');  
matlab2tikz('results/GaAs_p_doped_N_A_=1e21_chemical_potential_vs_temperature.tex','showInfo', false);
    
        
% intrinsic and main chrge carrier Si, Ge , GaAs%   
temperaturenew=linspace(10,1200,100);
[Si_p_dopednew] = simulate_semiconductor('Si', 0.1, 1e21,'p-doped', temperaturenew, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
[Ge_p_dopednew] = simulate_semiconductor('Ge', 0.1, 1e21,'p-doped', temperaturenew, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
[GaAs_p_dopednew] = simulate_semiconductor('GaAs', 0.1, 1e21,'p-doped', temperaturenew, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);

clf(figure(6))
figure(6)
hold on
set(gcf,'color','w');

plot(Si_p_dopednew.temperature,  Si_p_dopednew.n_i/Si_p_dopednew.dopant_density,'.','LineWidth',2,'Color',...
 [1 0 0],'DisplayName','Si_p_dopednew');
 plot(Ge_p_dopednew.temperature,  Ge_p_dopednew.n_i/Ge_p_dopednew.dopant_density,'.','LineWidth',2,'Color',...
 'g','DisplayName','Ge_p_dopednew');
 plot(GaAs_p_dopednew.temperature,  GaAs_p_dopednew.n_i/GaAs_p_dopednew.dopant_density,'.','LineWidth',2,'Color',...
 'b','DisplayName','GaAs_p_dopednew');
plot(Si_p_dopednew.temperature, Si_p_dopednew.main_charge_carrier_number/Si_p_dopednew.dopant_density,'LineWidth',2,'Color',...
 [1 0 0],'DisplayName','Si_p_dopednew');
 plot(Ge_p_dopednew.temperature, Ge_p_dopednew.main_charge_carrier_number/Ge_p_dopednew.dopant_density,'LineWidth',2,'Color',...
 'g','DisplayName','Ge_p_dopednew');
 plot(GaAs_p_dopednew.temperature, GaAs_p_dopednew.main_charge_carrier_number/GaAs_p_dopednew.dopant_density,'LineWidth',2,'Color',...
 'b','DisplayName','GaAs_p_dopednew');
     
ylim([0 2])
% title({'electron density vs temperature',' ',...
%           'in Si-p-doped, Ge-p-doped and GaAs-p-doped  at N_A = 10^{21} m^3'});
ax = gca;
ax.FontSize = 11;
grid on
legend('{(n_i/N_A)}_{Si}','{(n_i/N_A)}_{Ge}','{(n_i/N_A)}_{GaAs}','{(p/N_A)}_{Si}', ...
       '{(p/N_A)}_{Ge}','{(p/N_A)}_{GaAs}', 'Location' ,'northeastoutside');

  
xlabel('temperature / K');
ylabel('hole density / N_A');
matlab2tikz('results/All_p_doped_N_A_=1e21_ni_n_vs_NA.tex','showInfo', false);
    

% ionized dopants Si, Ge ,GaAs   
clf(figure(7))
figure(7)
% correct for evalation errors
set(gcf,'color','w');
Si_p_doped.ionized_dopants(find(Si_p_doped.ionized_dopants == 1)) = 0;
Ge_p_doped.ionized_dopants(find(Ge_p_doped.ionized_dopants == 1)) = 0;
GaAs_p_doped.ionized_dopants(find(GaAs_p_doped.ionized_dopants == 1)) = 0;
hold on
    
plot(Si_p_doped.temperature, Si_p_doped.ionized_dopants,'LineWidth',2,'Color',[1 0 0],...
     'DisplayName','Si');
plot(Ge_p_doped.temperature, Ge_p_doped.ionized_dopants,'LineWidth',2,'Color','g',...
     'DisplayName','Ge');
plot(GaAs_p_doped.temperature, GaAs_p_doped.ionized_dopants,'LineWidth',2,'Color','b',...
     'DisplayName','GaAs');
legend('{(N_A^-/N_A)}_{Si}','{(N_A^-/N_A)}_{Ge}','{(N_A^-/N_A)}_{GaAs}', 'Location' ,'northeastoutside');
    
% title({'number of ionized dopants vs temperature',' ',...
%            'in Si-p-doped, Ge-p-doped and GaAs-p-doped  at N_A = 10^{21} m^3'}); 

ax = gca;
ax.FontSize = 11;
grid on
    
ylim([0 1.1]);
xlabel('temperature / K');
ylabel('density of ionized dopants/ N_A');
matlab2tikz('results/All_p_doped_nof_ionized_dopants_vs_temperature.tex','showInfo', false);
    
   
%%%                         TASK 1(d)                                      %%%     

bandgap_multiplier_vec=linspace(0.1,2,5);
legendInfo = cell(length(bandgap_multiplier_vec),1);
bandgap_container = cell(1,length(bandgap_multiplier_vec));
for i = 1 : length(bandgap_multiplier_vec)
    bandgap_container{1,i} = simulate_semiconductor('Si', 0.1, 10^21,'p-doped', temperature, 0, 0, 0, bandgap_multiplier_vec(i), e_m_n_multiplier, e_m_p_multiplier);
end
bandgap_potential_vec= zeros([ length(bandgap_multiplier_vec) length(temperature)]);
for i = 1 : length(bandgap_multiplier_vec)
    bandgap_potential_vec(i,:) = bandgap_container{i}.chemical_potential;
end

figure(11)
set(gcf,'color','w');
hold on
for i=1: length(bandgap_multiplier_vec)
    plot(temperature,bandgap_potential_vec(i,:),'LineWidth',2)
    legendInfo{i} = ['E_g =' num2str(bandgap_multiplier_vec(i)) ' * E_{g_{Si}}' ];
end
hold off
grid on
xlabel('temperature / K');
ylabel('\mu/ eV');
legend(legendInfo,'Location','Northwest');
matlab2tikz('results/Si_p_doped_mu_vs_temperature_bangap_modified.tex','showInfo', false);


e_m_p_multiplier_vec=linspace(0.1,2,5);
legendInfo = cell(length(e_m_p_multiplier_vec),1);
mass_container = cell(1,length(e_m_p_multiplier_vec));
for i = 1 : length(e_m_p_multiplier_vec)
    mass_container{1,i} = simulate_semiconductor('Si', 0.1, 10^21,'p-doped', temperature, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier_vec(i));
end
mass_potential_vec= zeros([ length(e_m_p_multiplier_vec) length(temperature)]);
for i = 1 : length(e_m_p_multiplier_vec)
    mass_potential_vec(i,:) = mass_container{i}.chemical_potential;
end

figure(12)
set(gcf,'color','w');
hold on
for i=1: length(e_m_p_multiplier_vec)
    plot(temperature,mass_potential_vec(i,:),'LineWidth',2)
    legendInfo{i} = ['e\_m\_p = ' num2str(e_m_p_multiplier_vec(i)) ' * e\_m\_n_{Si}' ];
end
hold off
grid on
xlabel('temperature / K');
ylabel('\mu/ eV');
legend(legendInfo,'Location','Northwest');
matlab2tikz('results/Si_p_doped_mu_vs_temperature_e_m_p_modified.tex','showInfo', false);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         TASK 2                                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------- common parameters for task 2 -----------------------
temperature = 300;
ND_vector = logspace(13,25,100);
donor_state_energy = Si_p_doped.E_C - 0.1;

%############################## a #########################################

%-------------------------- simulation loop -> task a) --------------------
data_container_a = cell(1,length(ND_vector));
for i = 1 : length(ND_vector)
    data_container_a{1,i} = simulate_semiconductor('Si', donor_state_energy, ND_vector(i),'n-doped', temperature, 0, 0, 0, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
end

%---------------------- extract needed data from container -> task a) -----
chemical_potential_vector_a = zeros([1 length(ND_vector)]);
chemical_potential_vector_intrinsic = zeros([1 length(ND_vector)]);
for i = 1 : length(ND_vector)
    chemical_potential_vector_a(i) = data_container_a{i}.chemical_potential;
    chemical_potential_vector_intrinsic(i) = data_container_a{i}.chemical_potential_i;
end

%############################## b #########################################

%-------------------------- simulation loop -> task b) --------------------
trap_state_energy = (data_container_a{1}.E_C / 2 - 0.1);
trap_state_density = 1e22;
trap_state_type = 'n-charged';
ND_vector = logspace(13,25,100);
data_container_b = cell(1,length(ND_vector));
for i = 1 : length(ND_vector)
    data_container_b{1,i} = simulate_semiconductor('Si', donor_state_energy, ND_vector(i),'n-doped', temperature, trap_state_energy, trap_state_density, trap_state_type, bandgap_multiplier, e_m_n_multiplier, e_m_p_multiplier);
end

%---------------------- extract needed data from container -> task b) -----
chemical_potential_vector_b = zeros([1 length(ND_vector)]);
for i = 1 : length(ND_vector)
    chemical_potential_vector_b(i) = data_container_b{i}.chemical_potential;
end

%###################### Results of Task 2 a) and b) ####################### 

%------------------------ plot combined results -> task a) and b) ---------
figure();
set(gcf,'color','w');
semilogx(ND_vector, chemical_potential_vector_a,'g-.','LineWidth',2);
hold on;
semilogx(ND_vector, chemical_potential_vector_b,'cyan','LineWidth',2);
semilogx(ND_vector,Si_p_doped.E_V * ones(size(ND_vector)),... 
'LineWidth',2,'Color',[0 0 1],'DisplayName','Si E_V');
semilogx(ND_vector, chemical_potential_vector_intrinsic,'--','LineWidth',2,...
         'Color',[1 0 0],'DisplayName','Si E_F_intrinsic');
semilogx(ND_vector,Si_p_doped.E_C * ones(size(ND_vector)),... 
'LineWidth',2,'Color','k','DisplayName','Si E_C');
hold off;
ax = gca;
ax.FontSize = 11;
grid on
xlabel('N_D / m^{-3}');
ylabel('\mu / ev');
xlim([ND_vector(1), ND_vector(end)]);
legend('\mu','\mu_t','E_V','\mu_i','E_C', 'Location' ,'northeastoutside');
matlab2tikz('results/Si_with_trapped_states.tex','showInfo', false);

rmpath('matlab2tikz/src')
