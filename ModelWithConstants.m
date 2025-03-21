%%
clear;
close all;
clc;
opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

%% Defining global variables
global r_on r V_1 V_3 V_m K_10 K_12 K_21 K_m K_34 K_43 K_30;

%% Constants
% From D'Argenio:
% V_150, V_3=15, V_m=1,
% K_10=0.04, K_12=0.2, K_21=0.1, K_m=10, K_34=0.9, K_43=0.4, K_30=2

r_on = 1000;    % Input when the dose is being given
r = r_on;       % Input
V_1 = 50;       % Volume of the bloodstream
V_3 = 15;       % Volume of the target cell
V_m = 1;        % Volume of stuff moving from bloodstream to target cell
% Constants between each compartment
K_10 = 0.04;
K_12 = 0.2;
K_21 = 0.1;
K_m = 10;
K_34 = 0.9;
K_43 = 0.4;
K_30 = 2;

%% Initial conditions
Tf = 40;
% Amount in each compartment at time 0
x1_0 = 0; 
x2_0 = 0;
x3_0 = 0;
x4_0 = 0;
X0 = [x1_0,x2_0,x3_0,x4_0];

%% Call The function for ode
[T,Xf] = ode45(@PKPD, [0,Tf], X0, opts);
y1 = Xf(:,1)/V_1;
y2 = Xf(:,3)/V_3;


%% ODE function
% For indices: X[ 1,  2,  3,  4]
% The values:   [x1, x2, x3, x4]
function dX = PKPD(t,X)
    global r_on r V_1 V_3 V_m K_10 K_12 K_21 K_m K_34 K_43 K_30;
    % Dosage calculation
    % if(mod(t,10) < 0.1)
    %     r=r_on;
    % end
    if(mod(t,10) >= 0.1)
        r = 0;
    end
    % Calculate the dX/dt
    dX = zeros(length(X),1);
    dX(1) = -((V_1*V_m)/(V_1*K_m+X(1))+K_10+K_12)*X(1) + K_21*X(2)+r;
    dX(2) = K_12*X(1) - K_21*X(2);
    dX(3) = ((V_1*V_m)/(V_1*K_m+X(1)))*X(1) - (K_30+K_34)*X(3) + K_43*X(4);
    dX(4) = K_34*X(3) - K_43*X(4);
end

figure(1);
plot(T,Xf(:,1),'-r','LineWidth',2); hold on;
plot(T,Xf(:,2),'-g','LineWidth',2);
plot(T,Xf(:,3),'-b','LineWidth',2);
plot(T,Xf(:,4),'-k','LineWidth',2);
legend('Blood Plasma', 'Erythrocytes', 'Target Cell', 'Nucleus');
title('PK/PD Model of Doxorubicin Mass in the body over time');
xlabel('Time t (hours)');
ylabel('Amount of Doxorubicin in container (mg)');
set(gca,'fontsize',20);
set(gcf,'color','w');

figure(2);
plot(T,y1,'-m','LineWidth',2); hold on;
plot(T,y2,'-c','LineWidth',2);
legend('Blood Plasma','Target Cell');
title('PK/PD Model of Doxorubicin Concentration in the body over time');
xlabel('Time t (hours)');
ylabel('Concentration of Doxorubicin in container (mg/mL)');
set(gca,'fontsize',20);
set(gcf,'color','w');