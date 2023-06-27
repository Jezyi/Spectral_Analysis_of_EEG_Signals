%% compute empirical confidence intervals of AR estimator of the PSD of a simulated process
clear; close all; clc

%%% Parameters
N=300; % length of simulated time series (increase to see well autocorrelation)
numsimu=200;%quante volte riproduco il processo per stimare bias e varianza
fs=1;

SelCrit='BIC'; % model order selection criterion: AIC or BIC
pmax=10; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)


%% generate simu and compute PSD
par.poles{1}=([0.65 0; 0.8 0.1; 0.92 0.25]); % Oscillations RR
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes

 %%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(1,par); % parameters
%%% Theoretical (true) PSD
out_t = PSDvar(Am,Su,nfft,fs);
Px=squeeze(out_t.P);

%%% estimation
for n=1:numsimu
    clc; disp(['simulation ' int2str(n) ' of ' int2str(numsimu)])
    
    U=sqrt(Su)*randn(N,1);
    X=generateVAR(Am,U);

    %%% Estimated PSD
    % model order selection
    ret=VARorder(X,pmax);
    % [pottaic,pottmdl,aic,bic] = mos_idMVAR(S',pmax,0); [ret.aic aic ret.bic bic]
    switch SelCrit
        case 'AIC'
            ep=ret.pottaic;
        case 'BIC'
            ep=ret.pottbic;
    end
    % AR model identification
    jv=1; % index of predicted series
    iv=1; % indexes of predictors
    iv_lags=(1:ep); % lags of predictors
    out_e=LinReg(X,jv,iv,iv_lags);
    eAm=out_e.eA;
    eSu=out_e.es2u;
    % spectrum estimation
    out_e = PSDvar(eAm,eSu,nfft,fs);
    ePx(:,n)=squeeze(abs(out_e.P));
    f=out_e.f';
end

ePxm=mean(ePx,2); % mean PSD
ePxsd=std(ePx,0,2); % std. dev. PSD


%% plots and disps
colm=[86 86 158]/255; % color of mean value
colsh=[189 185 219]/255; %color of shades

figure(1);
hold on;
plot(f,Px,'k--','LineWidth',1.5);
plot(f,ePxm,'Color',colm,'LineWidth',1.5);
h1=area(f,ePxm+ePxsd);
h2=area(f,ePxm-ePxsd);
plot(f,Px,'k--','LineWidth',1.5);
plot(f,ePxm,'Color',colm,'LineWidth',1.5);
legend('true PSD','est.PSD, mean','est.PSD, +-SD')
set(h2,'FaceColor','w');
set(h1,'FaceColor',colsh);
set(h1,'EdgeColor',colsh); set(h2,'EdgeColor','w');
xlim([0 fs/2])
xlabel('f')
ylabel('P_x(f)')
