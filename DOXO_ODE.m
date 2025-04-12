%% A PK/PD model for Doxorubicin based on the paper by Daniele Andrean
clear;
close all;
clc;

%% Initial Conditions
Ihrs = 3;           %h
Itotl = 200;        %nM
Ion = Itotl/Ihrs;   %nM
I0 = 0;     %nM
X_E0 = I0;  %nM
X_F0 = 0;   %nM
X_B0 = 0;   %nM
N0 = 1e5;   %#
T0 = 0;     %h
Tf = 1000;  %h
% dt = 1e-4;
% nsteps = 1000/dt;
% opts = odeset(MinStep=dt, MaxStep=dt);
% opts = odeset(InitialStep=dt);

%% ODE function
X0 = [X_E0, X_F0, X_B0, N0];
[T,Xf] = ode45(@PKPD, [T0,Tf], X0);

function dX = PKPD(t,X)
    % Coefficients
    Ihrs = 3;           %h
    Itotl = 200;        %nM
    Ion = Itotl/Ihrs;   %nM
    % Constants
    V_max = 1.65e4; %nM/h
    k_th = 464;     %nM
    k_FE = 5.63e-4; %1/h
    k_BF = 1.22;    %1/h
    V_E = 100;      %µL
    V_I = 0.07;     %µL
    hl = 30;        %h (half-life)
    % Distributions (Normal????)
    k_p = normrnd(0.0198,0.002);        %1/h
    theta = normrnd(262209,13686);      %Dimensionless
    K_dmax = normrnd(0.0435,0.0025);    %1/h
    X_BHS = normrnd(47.5,9.83);         %nM
    gamma = normrnd(0.0044,0.0004);     %1/h
    if (t<50 || t>=53)
        I = 0;
    else
        I = Ion;
    end
    % Variables
    X_E = X(1);
    X_F = X(2);
    X_B = X(3);
    N = X(4);
    dX = zeros(length(X),1);
    % Calculate K_EF
    k_EF = V_max*X_B^1.31/(k_th^2.31+X_B^2.31);
    % Calculate k_d
    k_d = K_dmax*X_B/(X_BHS+X_B);
    % Calculate the concentrations of each compartment
    dX(1) = k_EF*V_I/V_E*X_F - k_FE*X_E + I - log(2)/hl*X_E; % - log(2)/hl*X_E added for half-life
    dX(2) = k_FE*V_E/V_I*X_E - k_EF*X_F - k_BF*X_F;
    dX(3) = k_BF*X_F - gamma*X_B - k_d*X_B; % - k_d*X_B added for cell death
    dX(4) = k_p*N*(1-N/theta) - k_d*N;
end

figure(1);
plot(T,Xf(:,1),'-r','LineWidth',2); hold on;
plot(T,Xf(:,2),'-g','LineWidth',2);
plot(T,Xf(:,3),'-b','LineWidth',2); hold off;
legend('X_E', 'X_F', 'X_B');
title('Pharmacokinetic Model of Doxorubicin');
xlabel('Time (hours)');
ylabel('Amount of Doxorubicin in container (mg)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(2);
plot(T,Xf(:,2)+Xf(:,3), '-m','linewidth',2); hold off;
legend('X_I');
title('Pharmacokinetic Model of Doxorubicin Concentration in the body over time');
xlabel('Time (hours)');
ylabel('Concentration of Doxorubicin (nM)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(3);
plot(T,Xf(:,4),'-k','LineWidth',2); hold off;
legend('N');
title('Pharmacodynamic Model of Doxorubicin');
xlabel('Time (hours)');
ylabel('Number of Cells');
set(gca,'fontsize',20);
set(gcf,'color','w');