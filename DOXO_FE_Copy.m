%% A PK/PD model for Doxorubicin based on the paper by Daniele Andrean
% Modeling a bolus injection of doxorubicin in the body and its effect
% on breast cancer (from an in vitro/in silico study)
clear;
close all;
clc;

%% Initial conditions
Ihrs = 3;           %h
Itot = 200;         %nM
Ion = Itot/Ihrs;    %nM
I0 = 0;     %nM
X_E0 = I0;  %nM
X_F0 = 0;   %nM
X_B0 = 0;   %nM
X_I0 = X_F0 + X_B0;
N0 = 1e5;   %#

%% Setting up the functions
T0 = 0;     %h
Tf = 1000;  %h
dt = 1e-4;
nsteps = Tf/dt;
T = linspace(T0,Tf,nsteps);
I = I0;
X_E = zeros(1,nsteps); X_E(1) = X_E0;
X_F = zeros(1,nsteps); X_F(1) = X_F0;
X_B = zeros(1,nsteps); X_B(1) = X_B0;
X_I = zeros(1,nsteps); X_I(1) = X_I0;
N   = zeros(1,nsteps); N(1)   = N0;

%% Coefficient Values
% Constants
dhl = 5/60;     %h (distribution half-life)
ehl = 48;       %h (elimination half-life)
hl = ehl;
q = 2.31;       %Dimensionless constant
V_max = 1.65e4; %nM/h
k_th = 464;     %nM
k_FE = 5.63e-4; %1/h
k_BF = 1.22;    %1/h
V_E = 100;      %µL
V_I = 0.07;     %µL
% Distributions
k_p = normrnd(0.0198,0.002);        %1/h
theta = normrnd(262209,13686);      %Dimensionless
K_dmax = normrnd(0.0435,0.0025);    %1/h
X_BHS = normrnd(74.5,9.83);         %nM
gamma = normrnd(0.0044,0.004);     %1/h

%% Forward Euler
for i=1:nsteps-1
    % Calculate input
    t = i*dt;
    if ((t>50 && t<53))
        I = Ion; hl = dhl;
    elseif ((t==53))
        I=0; hl = ehl;
        X_E(i) = 0;
    end

    % Calculate K_EF
    k_EF = V_max*((X_B(i))^1.31)/(k_th^2.31+(X_B(i))^2.31); % Exponents: q-1, q, q
    % Calculate k_d
    k_d = K_dmax*X_B(i)/(X_BHS+X_B(i));
    % Calculate the concentrations of each compartment
    X_E(i+1) = X_E(i) + dt*(k_EF*(V_I/V_E)*X_F(i) - k_FE*X_E(i) + I - log(2)/hl*X_E(i)); % - log(2)/hl*X_E(i) added for half-life
    X_F(i+1) = X_F(i) + dt*(k_FE*(V_E/V_I)*X_E(i) - k_EF*X_F(i) - k_BF*X_F(i));
    X_B(i+1) = X_B(i) + dt*(k_BF*X_F(i) - gamma*X_B(i) - k_d*X_B(i)); % - k_d*X_B(i) added for cell death
    X_I(i+1) = X_F(i+1) + X_B(i+1);
    N(i+1)   = N(i)   + dt*(k_p*N(i)*(1-(N(i)/theta)) - k_d*N(i));
end


%% Plots
% figure(1);
% plot(T,X_E, '-g','linewidth',2); hold off;
% legend('X_E');
% title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
% xlabel('Time (hours)');
% ylabel('Concentration of Doxorubicin (nM)');
% set(gca,'fontsize',20);
% set(gcf,'color','w');
% 
% figure(2);
% plot(T,X_F, '-g','linewidth',2); hold on;
% plot(T,X_B, '-b','linewidth',2); hold off;
% legend('X_F','X_B');
% title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
% xlabel('Time (hours)');
% ylabel('Concentration of Doxorubicin (nM)');
% set(gca,'fontsize',20);
% set(gcf,'color','w');

figure(3);
plot(T,X_I, '-r','linewidth',2); hold off;
legend('X_I');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time (hours)');
ylabel('DOXO (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(4);
plot(T,N, '-m','linewidth',2); hold off;
legend('N');
title('Pharmacodynamic Model of Doxorubicin in the body over time');
xlabel('Time (hours)');
ylabel('Number of Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');
axis([0,600,0,2e5]);