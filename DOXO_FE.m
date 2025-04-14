%% A PK/PD model for Doxorubicin based on the paper by Daniele Andrean
% Modeling a bolus injection of doxorubicin in vitro and its effect
% on breast cancer cells (MM1R) in vitro
clear;
close all;
clc;

for Itot = [50,100,200] %nM
%% Initial conditions
% Duration of dose Ihrs, Amount per hour Ion
Ihrs = 3;           %h
Ion = Itot;    %nM/h
% Initial conditions
I0 = 0;     %nM
X_E0 = I0;  %nM
X_F0 = 0;   %nM
X_B0 = 0;   %nM
X_I0 = X_F0 + X_B0;
N0 = 1e5;   %#

%% Setting up the functions
% Time start, end, step, number of steps
T0 = 0;     %h
Tf = 1000;  %h
dt = 1e-4;  %h
nsteps = Tf/dt;
% Set up I and the Vectors
I = I0;
T = linspace(T0,Tf,nsteps);
X_E = zeros(1,nsteps); X_E(1) = X_E0;
X_F = zeros(1,nsteps); X_F(1) = X_F0;
X_B = zeros(1,nsteps); X_B(1) = X_B0;
X_I = zeros(1,nsteps); X_I(1) = X_I0;
N   = zeros(1,nsteps); N(1)   = N0;

%% Coefficient Values
% Constants
V_max = 1.65e4; %nM/h
k_th = 464;     %nM
k_FE = 5.63e-4; %1/h
k_BF = 1.22;    %1/h
V_E = 100;      %µL
V_I = 0.07;     %µL

% Distributions
% k_p = normrnd(0.0198,0.002);        %1/h
% theta = normrnd(262209,13686);      %Dimensionless
% K_dmax = normrnd(0.0435,0.0025);    %1/h
% X_BHS = normrnd(74.5,9.83);         %nM
% gamma = normrnd(0.0044,0.004);      %1/h

k_p = 0.0198;       %1/h
theta = 262209;     %Dimensionless
K_dmax = 0.0435;    %1/h
X_BHS = 47.5;       %nM
gamma = 0.0044;     %1/h

%% Forward Euler
for i=1:nsteps-1
    % Calculate input
    t = i*dt;
    if (Itot == 200)
        if (t==50)
            X_E(i) = Ion;
        elseif ((t==53))
            X_E(i) = 0;
        end
    elseif (Itot == 100)
        if ((t==50) || (t==350))
            X_E(i) = Ion;
        elseif ((t==53) || (t==353))
            X_E(i) = 0;
        end
    elseif (Itot == 50)
        if ((t==50) || (t==150) || (t==300) || (t==450))
            X_E(i) = Ion;
        elseif ((t==53) || (t==153) || (t==303) || (t==453))
            X_E(i) = 0;
        end
    end

    % Calculate K_EF
    k_EF = V_max*((X_B(i))^1.31)/(k_th^2.31+(X_B(i))^2.31); %Exponents: q-1, q, q
    % Calculate k_d
    k_d = K_dmax*X_B(i)/(X_BHS+X_B(i));
    % Calculate the concentrations of each compartment
    X_E(i+1) = X_E(i) + dt*(k_EF*(V_I/V_E)*X_F(i) - k_FE*X_E(i) + I);
    X_F(i+1) = X_F(i) + dt*(k_FE*(V_E/V_I)*X_E(i) - k_EF*X_F(i) - k_BF*X_F(i));
    X_B(i+1) = X_B(i) + dt*(k_BF*X_F(i) - gamma*X_B(i));
    X_I(i+1) = X_F(i+1) + X_B(i+1);
    N(i+1)   = N(i)   + dt*(k_p*N(i)*(1-(N(i)/theta)) - k_d*N(i));
end
% Save the Vectors for Each Itot
if (Itot == 200)
    X_F_200 = X_F;
    X_B_200 = X_B;
    X_I_200 = X_I;
    N_200 = N;
elseif (Itot == 100)
    X_F_100 = X_F;
    X_B_100 = X_B;
    X_I_100 = X_I;
    N_100 = N;
elseif (Itot == 50)
    X_F_50 = X_F;
    X_B_50 = X_B;
    X_I_50 = X_I;
    N_50 = N;
end

%% Plots
% % X_E graph for current Itot
% figure(1);
% plot(T,X_E, '-g','linewidth',2); hold off;
% legend('X_E');
% title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
% xlabel('Time (hours)');
% ylabel('Concentration of Doxorubicin (nM)');
% set(gca,'fontsize',20);
% set(gcf,'color','w');
% 
% % X_F, X_B graph for current Itot
% figure(2);
% plot(T,X_F, '-g','linewidth',2); hold on;
% plot(T,X_B, '-b','linewidth',2); hold off;
% legend('X_F','X_B');
% title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
% xlabel('Time (hours)');
% ylabel('Concentration of Doxorubicin (nM)');
% set(gca,'fontsize',20);
% set(gcf,'color','w');
% 
% % X_I graph for current Itot
% figure(Itot+1);
% plot(T,X_I, '-r','linewidth',2); hold off;
% legend('X_I');
% title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
% xlabel('Time (hours)');
% ylabel('DOXO (nM)');
% set(gca,'fontsize',20);
% set(gcf,'color','w');
% 
% % N graph for current Itot
% figure(Itot+2);
% plot(T,N, '-m','linewidth',2); hold off;
% legend('N');
% title('Pharmacodynamic Model of Doxorubicin in the body over time');
% xlabel('Time (hours)');
% ylabel('Number of Cells');
% set(gca,'fontsize',20);
% set(gcf,'color','w');
% axis([0,600,0,2e5]);
end

%% Final Plots
% X_I graph for all Itot
figure(5);
plot(T,X_I_200, '-r','linewidth',2); hold on;
plot(T,X_I_100, '-g','linewidth',2); hold on;
plot(T,X_I_50, '-b','linewidth',2); hold off;
legend('200 nM', '100 nM', '50 nM');
title('Pharmacokinetic Model of Doxorubicin Concentration');
xlabel('Time (hours)');
ylabel('Concentration of Doxorubicin (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% First bolus zoom
figure(6);
plot(T,X_F_200, '-', 'color',[1 0 0 0.1], 'linewidth',2); hold on;
plot(T,X_F_100, '-', 'color',[0 1 0 0.1], 'linewidth',2); hold on;
plot(T,X_F_50, '-', 'color',[0 0 1 0.1], 'linewidth',2); hold on;
plot(T,X_B_200, '-', 'color',[1 0 0 0.3], 'linewidth',2); hold on;
plot(T,X_B_100, '-', 'color',[0 1 0 0.3], 'linewidth',2); hold on;
plot(T,X_B_50, '-', 'color',[0 0 1 0.3], 'linewidth',2); hold on;
plot(T,X_I_200, '-r','linewidth',2); hold on;
plot(T,X_I_100, '-g','linewidth',2); hold on;
plot(T,X_I_50, '-b','linewidth',2); hold off;
legend('X_F', '', '', 'X_B', '', '', 'X_I');
title('Pharmacokinetic Model of Doxorubicin Concentration');
xlabel('Time (hours)');
ylabel('Concentration of Doxorubicin (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');
axis([50,56,0,160]);

% X_B
figure(7);
plot(T,X_B_200, '-r','linewidth',2); hold on;
plot(T,X_B_100, '-g','linewidth',2); hold on;
plot(T,X_B_50, '-b','linewidth',2); hold off;
legend('200 nM', '100 nM', '50 nM');
title('Pharmacokinetic Model of Doxorubicin Concentration');
xlabel('Time (hours)');
ylabel('Concentration of Doxorubicin (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% N graph for all Itot
figure(8);
plot(T,N_200, '-r','linewidth',2); hold on;
plot(T,N_100, '-g','linewidth',2); hold on;
plot(T,N_50, '-b','linewidth',2); hold off;
legend('200 nM', '100 nM', '50 nM');
title('Pharmacodynamic Model of Breast Cancer Cells');
xlabel('Time (hours)');
ylabel('Number of Cells (#)');
set(gca,'fontsize',20);
set(gcf,'color','w');
axis([0,600,0,2e5]);