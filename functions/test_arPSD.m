%% test of AR PSD estimators on simulated RR time series
clear; close all; clc;

%%%% parameters
SelCrit='BIC'; % model order selection criterion: AIC or BIC
pmax=10; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)

%% simulation of RR process (time series x)
N=300; % length of simulated time series (increase to see well autocorrelation)
p=5; % model order
fs=1; % sampling frequency corresponding to mean RR=1 sec

par.poles{1}=([0.65 0; 0.8 0.1; 0.92 0.25]); % Oscillations RR
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes
%%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(1,par); % parameters

% Realization of the simulation
U=sqrt(Su)*randn(N,1);
x=generateVAR(Am,U);

%% Power Spectral Density
%%% Theoretical (true) PSD
out_t = PSDvar(Am,Su,nfft,fs);
Px=squeeze(out_t.P);

%%% Estimated PSD
% model order selection
ret=VARorder(x,pmax);
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
out_e=LinReg(x,jv,iv,iv_lags);
eAm=out_e.eA;
eSu=out_e.es2u;

% spectrum estimation
out_e = PSDvar(eAm,eSu,nfft,fs);
ePx=squeeze(abs(out_e.P));
f=out_e.f';

% to get the integral, each PSD point is multiplied for 2*pi/nfft; the integral is then divided by 2*pi to get the variance (see slide 20)
Ptot=sum(ePx)/nfft;

%% plots
disp(['model order, true:p=' int2str(p) '; estimated:p=' int2str(ep)])
disp('AR parameters [true estimated]:')
disp([[Am; Su] [eAm; eSu]])

figure(1);
subplot(2,1,1); plot(x); title(['time series, variance:' num2str(var(x))])
xlabel('n'); ylabel('x_n')
subplot(2,1,2);
plot(f,Px,'k--','LineWidth',1.2);
hold on; 
plot(f,ePx,'LineWidth',1.5);
xlim([0 fs/2])
title(['Power spectral density, Ptot:' num2str(Ptot)])
xlabel('f'); ylabel('P_X(f)')
legend('true PSD','estimated PSD')



