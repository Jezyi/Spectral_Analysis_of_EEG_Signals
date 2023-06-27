%% calculation of the PSD of a representative EEG
clear; close all; clc;

subjfold='01'; % subject (name of folder) - from 01 to 21
condition='aperti'; % 'aperti' or 'chiusi'
electrode=19;  % EEG to analyze (channel): from 1 to 19

falpha=[8 12]; % limits of alpha band (Hz)

%%%% PSD parameters
winname='blackman'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='unbiased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
M=120; % lag at which to truncate the correlation
nfft=1000; %number of points on frequency axis (total)

%% load EEG data and show one signal
directory = '..\exe_01_Database\D-EEG\data\';
load([directory subjfold '\' condition '\sottocamp_data_viever.mat']);
data=v_salvato.ser;
chan=v_salvato.nomec';
fs=v_salvato.fc;

x = data(:,11+electrode);
ch_name = chan(11+electrode);


%% spectral analysis
[Px,f,Pxtot,ftot,Omega]=PSDwc(x,M,typecorrest,winname,fs,nfft);
N=length(x);% number of time series points
Nf=length(Px); %number of frequency points

% total EEG power: var(x)=(2/fs)*(integral of Px from 0 to fs/2);
% (2/fs) * (sum(Px) * (fs/2)/Nf ) = sum(Px)/Nf
Ptot=sum(Px)/Nf;

% f:n=(fs/2):Nf
nalpha=round(falpha*Nf*2/fs); % range alpha in frequency bins
if nalpha(1)==0, nalpha(1)=1; end %adjust inferior bound if needed
if nalpha(2)==Nf, nalpha(2)=Nf; end %adjust superior bound if needed
Palpha=sum(Px(nalpha(1):nalpha(2)))/Nf;
Palpha_n=Palpha/Ptot;

%% plots
t=(1/fs:1/fs:N/fs)'; % time axis
colsh=[189 185 219]/255; %color of shades

figure(1)
subplot(2,1,1)
plot(t,x,'color',[0 0.4 0.7]);
xlim([t(1) t(N)])
xlabel('time [sec]');
% text(t(10),0.95*max(x),['variance = ' num2str(var(x))])
title(['EEG signal Electrode ', char(ch_name) ', ' condition])

subplot(2,1,2)
plot(f,Px,'color',[0 0.4 0.7],'LineWidth',1.5);
hold on
h1=area(f(nalpha(1):nalpha(2)),Px(nalpha(1):nalpha(2)));
set(h1,'FaceColor',colsh);
title(['PSD, Electrode ', char(ch_name) ', ' condition])
xlim([0 fs/2])
xlabel('frequency [Hz]');


disp(['variance = ' num2str(var(x))])
disp(['total power = ' num2str(Ptot)])
disp(['alpha power = ' num2str(Palpha)])
disp(['alpha power, normalized = ' num2str(100*Palpha_n) ' %'])


