%% compute empirical confidence intervals of WC estimator of the PSD of a simulated process
clear; close all; clc

%%% Parameters
N=300; % length of simulated time series (increase to see well autocorrelation)
numsimu=200;%quante volte riproduco il processo per stimare bias e varianza

winname='blackman'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('hamming', 'hanning',...)
typecorrest='unbiased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
M=40; % lag at which to truncate the correlation
nfft=1001; %number of points on frequency axis (total)
fs=1;


%% generate simu and compute PSD
par.poles{1}=([0.65 0; 0.8 0.1; 0.92 0.25]); % Oscillations RR
par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=1; %variance of innovation processes

 %%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(1,par); % parameters
for n=1:numsimu
    U=sqrt(Su)*randn(N,1);
    X=generateVAR(Am,U);

    [ePxn,f,ePx2n]=PSDwc(X,M,typecorrest,winname,fs,nfft);
    ePx(:,n)=ePxn;
%     ePx2(:,n)=ePx2n;
end

ePxm=mean(ePx')'; % mean PSD
ePxsd=std(ePx')'; % std. dev. PSD

% spectral analysis and true spectrum
% Apol=[1 -Am']; [Px,ft]=ar_spectrum(Apol,Su,ceil(nfft/2),fs);
out_t = PSDvar(Am,Su,ceil(nfft/2),fs); Px=squeeze(out_t.P); ft=out_t.f;

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
