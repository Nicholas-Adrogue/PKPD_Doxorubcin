%%
clear;
close all;
clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% Defining global variables
global tox int_hrs r_on r V_1 V_2 V_3 V_4 V_m CC CellRep K_10 K_12 K_21 K_m K_34 K_43 K_30;

%% Constants
% From D'Argenio:
% V_150, V_3=15, V_m=1,
% K_10=0.04, K_12=0.2, K_21=0.1, K_m=10, K_34=0.9, K_43=0.4, K_30=2

% Dosage Calculation
dose = 60;       % (mg/m^2)
if (rand()<0.5)  % Male or Female height (cm)
    Height = normrnd(178.4,7.59);   % Male
else
    Height = normrnd(164.7,7.07);   %Female
end
Weight = 70;     % (kg)
BSA = (Height*Weight/3600)^0.5; % Mosteller's Body Surface Area (m^2)
r_on = dose*BSA;    % Input when the dose is being given (mg)

% Constants
r = r_on;       % Input
V_1 = 5000;     % Volume of the bloodstream (mL)
V_2 = 90*1e-12; % Volume of the erythrocytes (mL)
V_3 = 2.5e-8;   % Volume of the target cell (mL)
V_4 = 0.1*V_3;  % Volume of the nucleus of the target cell (mL)
V_m = 24*1000/24; % Volume of stuff moving from bloodstream to target cell (mL/Hr)
CC = 500000;    % Carrying capacity of heart cells (#)
CellRep = 0.01*CC/8760; % Speed at which heart cells are replaced (#/hour)
tox = 20; % Toxicity coefficient

% Constants between each compartment
K_10 = .04;
K_12 = .2;
K_21 = .1;
K_m = 10;
K_34 = .09;
K_43 = .04;
K_30 = 2;

%% Initial conditions
interval = 21; % (Days) Time between doses
int_hrs = 21*24; % (Hours) Time between doses
NumDoses = 4; % (#) Number of doses total
Tf = int_hrs*NumDoses; % (Hours) Time total
% Amount in each compartment at time 0
x1_0 = 0;
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
c_0 = CC; % Cells in human heart
X0 = [x1_0,x2_0,x3_0,x4_0,c_0];

%% Call The function for ode
[T,Xf] = ode45(@PKPD, [0,Tf], X0, opts);
y1 = Xf(:,1)/V_1;
y2 = Xf(:,3)/V_3;


%% ODE function
% For indices: X[ 1,  2,  3,  4,  5]
% The values:   [x1, x2, x3, x4,  c]
function dX = PKPD(t,X)
    global tox int_hrs r_on r V_1 V_2 V_3 V_4 V_m CC CellRep K_10 K_12 K_21 K_m K_34 K_43 K_30;
    % Dosage calculation
    if(mod(t,int_hrs) < 0.5)
        r=r_on;
    end
    if(mod(t,int_hrs) >= 0.5)
        r = 0;
    end
    % Calculate the dX/dt
    dX = zeros(length(X),1);
    dX(1) = -((V_1*V_m)/(V_1*K_m+X(1))+K_10+K_12)*X(1) + K_21*X(2)+r;
    dX(2) = K_12*X(1) - K_21*X(2);
    dX(3) = ((V_1*V_m)/(V_1*K_m+X(1)))*X(1) - (K_30+K_34)*X(3) + K_43*X(4);
    dX(4) = K_34*X(3) - K_43*X(4);
    dX(5) = CellRep*X(5)*(1-X(5)/CC) - CellRep*X(5)*(1-CC/X(5)) - tox*X(4);
end

figure(1);
plot(T,Xf(:,1),'-r','LineWidth',2); hold on;
plot(T,Xf(:,2),'-g','LineWidth',2);
plot(T,Xf(:,3),'-b','LineWidth',2);
plot(T,Xf(:,4),'-k','LineWidth',2); hold off;
legend('Blood Plasma', 'Erythrocytes', 'Target Cell', 'Nucleus');
title('Pharmacokinetic Model of Doxorubicin Mass in the body over time');
xlabel('Time t (hours)');
ylabel('Amount of Doxorubicin in container (mg)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(2);
plot(T,y1,'-m','LineWidth',2); hold on;
plot(T,y2,'-c','LineWidth',2); hold off;
legend('Blood Plasma','Target Cell');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time t (hours)');
ylabel('Concentration of Doxorubicin in container (mg/mL)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(3);
plot(T,Xf(:,5),'-m','LineWidth',2);
title('Pharmacodynamic Model of Heart Cells in the body');
xlabel('Time t (hours)');
ylabel('Number of Heart Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');