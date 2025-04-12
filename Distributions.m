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
title('Max Dosage of DOXO');
xlabel('Dosage (mg/m^2/day)');
set(gca,'fontsize',20);
set(gcf,'color','w');

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
xlabel('Concentration (ng/mL)');
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
title('Average Clearance of DOXO');
xlabel('Clearance (ml/min/m^2)');
set(gca,'fontsize',20);
set(gcf,'color','w');