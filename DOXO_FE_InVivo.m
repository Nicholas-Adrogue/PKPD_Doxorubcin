%% A PK/PD model for Doxorubicin based on the paper by Ackland
% Modeling a continuous injection of doxorubicin in the body and its effect
% on breast cancer cells in vivo
clear;
close all;
clc;
PatData = readtable("synthetic_patient_data.csv");
ToxData = readtable("synthetic_toxicity_data.csv");

%% Dosage Distribution
R_I = PatData.MaximumInfusionRate_mg_m_3_day_;
Mn = mean(R_I);
Sd = std(R_I);
Sk = skewness(R_I);
Ku = kurtosis(R_I);

X = linspace(2,6,100);
DoseNorm = fitdist(R_I,'Normal');
DoseSkew = pearspdf(X,Mn,Sd,Sk,Ku);

figure(101);
histogram(R_I,'Normalization','pdf'); hold on;
plot(X,DoseSkew, 'LineWidth',2); hold off;
xlim([2,6]);

dose = pearsrnd(Mn,Sd,Sk,Ku,1); %mg/m^2/day

%% Max Concentration Distribution
MaxDox = ToxData.MaximumDoxorubicinConcentration_ng_ml_+ToxData.MaximumDoxorubicinolConcentration_ng_ml_;
MnDox = mean(MaxDox);
SdDox = std(MaxDox);
SkDox = skewness(MaxDox);
KuDox = kurtosis(MaxDox);

X = linspace(0,25,100);
DoxNorm = fitdist(MaxDox,'Normal');
DoxSkew = pearspdf(X,MnDox,SdDox,SkDox,KuDox);

figure(102);
histogram(MaxDox,'Normalization','pdf'); hold on;
plot(X,DoxSkew, 'LineWidth',2); hold off;
xlim([0,25]);
title('Max Concentration of DOXO and DOXol');
xlabel('Time (hours)');
ylabel('Number of Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');

Dox = pearsrnd(MnDox, SdDox, SkDox, KuDox, 1); %ng/mL

%% Clearance Distribution
ClrM = PatData.MeanClearance_ml_min_m_2_;
MnClr = mean(ClrM);
SdClr = std(ClrM);
SkClr = skewness(ClrM);
KuClr = kurtosis(ClrM);

X = linspace(200,1000,100);
ClrNorm = fitdist(ClrM,'Normal');
ClrSkew = pearspdf(X,MnClr,SdClr,SkClr,KuClr);

figure(103);
histogram(ClrM,'Normalization','pdf'); hold on;
plot(X,ClrSkew, 'LineWidth',2); hold off;
xlim([200,1000]);

Clr = pearsrnd(MnClr, SdClr, SkClr, KuClr, 1); %mL/min/m^2

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
I = dose*BSA/MM/24/V_E;         %Input (nM/h)
Cm = Dox*1000;                  %Max concentration (nM)
CT = Clr*BSA*60/1000;           %Clearance Time (h)

%% Initial conditions
% Initial conditions
X_E0 = 0;   %nM
X_F0 = 0;   %nM
X_B0 = 0;   %nM
X_I0 = X_F0 + X_B0;
N0 = 1e8;   %#

%% Setting up the functions
% Time start, end, step, number of steps
T0 = 0;     %h
Tf = 10000; %h
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
V_I = 2/1000*N0;%L

% Distributions
% k_p = normrnd(0.0198,0.002);        %1/h
% theta = normrnd(262209,13686);      %Dimensionless
% K_dmax = normrnd(0.0435,0.0025);    %1/h
% X_BHS = normrnd(74.5,9.83);         %nM
% gamma = normrnd(0.0044,0.004);      %1/h

k_p = 0.0198;       %1/h
theta = 8.4e6;      %#
K_dmax = 0.0435;    %1/h
X_BHS = 47.5;       %nM
gamma = 0.0044;     %1/h

%% Forward Euler
for i=1:nsteps-1
    % Calculate input
    t = i*dt;
    if (t>2000)
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

% Plots
% X_E graph for current Itot
figure(1);
plot(T,X_E, '-g','linewidth',2); hold off;
legend('X_E');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time (hours)');
ylabel('DOXO (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% X_F, X_B graph for current Itot
figure(2);
plot(T,X_F, '-g','linewidth',2); hold on;
plot(T,X_B, '-b','linewidth',2); hold off;
legend('X_F','X_B');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time (hours)');
ylabel('DOXO (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% X_I graph for current Itot
figure(3);
plot(T,X_I, '-r','linewidth',2); hold off;
legend('X_I');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time (hours)');
ylabel('DOXO (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

% N graph for current Itot
figure(4);
plot(T,N, '-m','linewidth',2); hold off;
legend('N');
title('Pharmacodynamic Model of Doxorubicin in the body over time');
xlabel('Time (hours)');
ylabel('Number of Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');