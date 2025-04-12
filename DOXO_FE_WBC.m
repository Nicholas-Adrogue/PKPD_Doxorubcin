%% A PK/PD model for Doxorubicin based on the paper by Ackland
% Modeling a continuous injection of doxorubicin in the body and its effect
% on white blood cells in vivo
clear;
close all;
clc;
PatData = readtable("synthetic_patient_data.csv");
ToxData = readtable("synthetic_toxicity_data.csv");

WBC0 = ToxData.InitialWBC_x10_3_l;  %#/e3/µL
dose = ToxData.MaximumInfusionRate_mg_m_3_day_; %mg/m^2/day
WBCn = ToxData.NadirWBC_x10_3_l;    %#/e3/µL
maxDox = ToxData.MaximumDoxorubicinConcentration_ng_ml_+ ...
        ToxData.MaximumDoxorubicinolConcentration_ng_ml_; %ng/mL
maxDoxCalc = zeros(length(maxDox),1);
WBCnCalc = zeros(length(WBCn),1);

for j=1:length(WBC0)
%% Dosage Calculation 
MM = 543.52*1000/1e9;   %Molar Mass of DOXO (mg/nmol)
V_E = 5;        %L
if (rand()<0.5) %Male or Female height (cm)
    Height = normrnd(178.4,7.59);   %Male
else
    Height = normrnd(164.7,7.07);   %Female
end
Weight = 70;     %(kg)
BSA = (Height*Weight/3600)^0.5; %Mosteller's Body Surface Area (m^2)
I = dose(j)*BSA/MM/24/V_E;      %Input (nM/h)
Cm = maxDox(j)*1000;            %Max concentration (nM)

%% Initial conditions
% Initial conditions
X_E0 = 0;   %nM
X_F0 = 0;   %nM
X_B0 = 0;   %nM
X_I0 = X_F0 + X_B0;
N0 = WBC0(j)*1e3*V_E*1e6;   %#

%% Setting up the functions
% Time start, end, step, number of steps
T0 = 0;     %h
Tf = 5000; %h
Weeks = 20; %Weeks
dt = 1e-3;  %h
nsteps = Tf/dt;
% Set up I and the Vectors
T = linspace(T0,Tf,nsteps);
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
V_I = 4.68e-13*N0; %L

% Distributions
% k_p = normrnd(0.0198,0.002);        %1/h
% theta = normrnd(262209,13686);      %Dimensionless
% K_dmax = normrnd(0.0435,0.0025);    %1/h
% X_BHS = normrnd(74.5,9.83);         %nM
% gamma = normrnd(0.0044,0.004);      %1/h

k_p = 0.0198;       %1/h
V_blood = 5;        %L
theta = 11000*1e6*V_blood;      %#
K_dmax = 0.0435;    %1/h
X_BHS = 47.5;       %nM
gamma = 0.0044;     %1/h

%% Forward Euler
for i=1:nsteps-1
    % Calculate input
    t = i*dt;
    if (t>Weeks*168)
        I=0;
    end
    % Calculate K_EF
    k_EF = V_max*((X_B(i))^1.31)/(k_th^2.31+(X_B(i))^2.31); %Exponents: q-1, q, q
    % Calculate k_d
    k_d = K_dmax*X_B(i)/(X_BHS+X_B(i));
    % Calculate the concentrations of each compartment
    X_E(i+1) = X_E(i) + dt*(k_EF*(V_I/V_E)*X_F(i) - k_FE*X_E(i) + I - (log(2)/hl)*X_E(i));
    X_F(i+1) = X_F(i) + dt*(k_FE*(V_E/V_I)*X_E(i) - k_EF*X_F(i) - k_BF*X_F(i));
    X_B(i+1) = X_B(i) + dt*(k_BF*X_F(i) - gamma*X_B(i));
    X_I(i+1) = X_F(i+1) + X_B(i+1);
    N(i+1)   = N(i)   + dt*(k_p*N(i)*(1-(N(i)/theta)) - k_d*N(i));
end
maxDoxCalc(j) = X_E(2000/dt); %nM
WBCnCalc(j) = min(N); %#
end

% Plots
% X_E graph for current Itot
figure(1);
plot(T,X_E, '-r','linewidth',2); hold off;
legend('X_E');
title('Pharmacokinetic Model of DOXO');
xlabel('Time (hours)');
ylabel('DOXO (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% N graph for current Itot
figure(5);
plot(T,N, '-r','linewidth',2); hold off;
legend('N');
title('Pharmacodynamic Model of Doxorubicin');
xlabel('Time (hours)');
ylabel('Number of Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(6);
scatter(maxDox,(maxDoxCalc/1e3),'filled'); hold on;
plot([1:25],[1:25],'--k', 'LineWidth',1.5); hold off;
title('Pharmacokinetic Model of Doxorubicin');
xlabel('Actual Maximum DOXO concentration (ng/mL)');
ylabel('Calculated Maximum DOXO concentration (ng/mL)');
set(gca,'fontsize',20);
set(gcf,'color','w');
legend('Max DOXO','y=x');
yaxis([0,15]);

figure(7);
scatter(WBCn,(WBCnCalc/1e3/V_E/1e6),'filled'); hold on;
plot([1:25],[1:25],'--k', 'LineWidth',1.5); hold off;
title('Pharmacokinetic Model of Doxorubicin');
xlabel('Actual Minimum WBC count (#/mL)');
ylabel('Calculated Minimum WBC Count (#/mL)');
set(gca,'fontsize',20);
set(gcf,'color','w');