%% test of PSD estimators on simulated RR time series
clear; close all; clc;

%%%% parameters
winname='blackman'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='biased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
M=30; % lag at which to truncate the correlation
nfft=1000; %number of points on frequency axis (total)

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
% N=size(x,1);

%% Power Spectral Density
% [Pxpos,fpos,Px,f,Omega,rx,rxw,wi,lags,nfft]=PSDwc(x); %default executes correlogram/periodogram
[Pxpos,fpos,Px,f,Omega,rx,rxw,wi,lags,nfft]=PSDwc(x,M,typecorrest,winname,fs,nfft);

% to get the integral, each PSD point is multiplied for 2*pi/nfft; the integral is then divided by 2*pi to get the variance (see slide 20)
Ptot=sum(Px)/nfft;



%% plots
figure(1);
subplot(2,1,1); plot(x); title(['time series, variance:' num2str(var(x))])
xlabel('n'); ylabel('x_n')
subplot(2,1,2);
plot(fpos,Pxpos/2,'LineWidth',1.5);
xlim([0 fs/2])
title(['Power spectral density, Ptot:' num2str(Ptot)])
xlabel('f'); ylabel('P_X(f)')

%%% spectral analysis and true spectrum
% Apol=[1 -Am']; [Pxt,ft]=ar_spectrum(Apol,Su,ceil(nfft/2),fs);
out_t = PSDvar(Am,Su,ceil(nfft/2),fs); Pxt=squeeze(out_t.P); ft=out_t.f;
hold on; plot(ft,Pxt/2,'k--','LineWidth',1.2);
legend('estimated PSD','true PSD')


figure(2);
subplot(1,2,1); 
plot(lags,rx); hold on; plot(lags,rxw);
xlim([min(lags) max(lags)])
xlabel('k'); ylabel('r_x(k)')
title(['est. correlation, window ' winname ', trunc. lag ' int2str(M)])
legend('original','windowed')
subplot(1,2,2); 
plot(Omega,Px,'k'); %hold on; plot(f,Px2)
xlabel('\Omega'); ylabel('P_X(\Omega)')
xlim([-pi pi])
title('estimated PSD')

figure(3)
subplot(1,2,1); 
lagsN=-(N-1):1:(N-1);
wiN=zeros(length(lagsN),1);
wiN(N-(length(lags)-1)/2:N+(length(lags)-1)/2)=wi;
plot(lagsN,wiN,'.-')
xlim([-(N-1) N-1])
xlabel('k'); ylabel('w(k)');
title('lag window')
subplot(1,2,2);
Wi=fft(wi,nfft); %DTFT of window
Wiflip=[flip(abs(Wi(1:floor(nfft/2)))); abs(Wi(1:ceil(nfft/2)))]; % to show the DTFT ifrom -pi to pi
plot(Omega,Wiflip)
xlim([-pi pi])
title('spectral window')
xlabel('f'); ylabel('W(f)');


