%% test of CSD estimators on simulated RR-SAP-RESP time series
clear; close all; clc;

fs=1; % sampling frequency corresponding to mean RR=1 sec
N=3000; % length of simulated time series

%%%% parameters
winname='blackman'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='biased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
M=30; % lag at which to truncate the correlation
nfft=1001; %number of points on frequency axis (total)


%% Simulation of cardiovascular and cardiorespiratory interactions
% p=5; % model order
Ms=2;
c=0.8; %coupling from respiration to RR
par.poles{1}=([0.9 0.25]); % Oscillations resp
par.poles{2}=([0.8 0.1]); % Oscillations RR
par.coup=[1 2 2 c]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[2 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(Ms,par); % parameters

% Realization of the simulation
U = mvnrnd(zeros(1,Ms),Su,N);
X=generateVAR(Am,U);


%% data
xy=X;

%%% Cross-PSD
[Cxy,Px,Py,Pxy,f]=COHwc(xy,M,typecorrest,winname,fs,nfft);

% to get the integral, each PSD point is multiplied for 2*pi/nfft; the integral is then divided by 2*pi to get the variance (see slide 20)
Ptotx=sum(Px)/length(Px);
Ptoty=sum(Py)/length(Py);

% cross-spectral analysis and true spectrum
out_t = PSDvar(Am,Su,ceil(nfft/2),fs);
MSCt=abs(out_t.COH).^2; %coherence
Cxyt=squeeze(MSCt(1,2,:)); ft=out_t.f;


%% plots
x=xy(:,1); y=xy(:,2);

figure(1);
subplot(2,4,[1 2]);
plot(x,'k'); title(['time series x, var:' num2str(var(x))])
xlabel('n'); ylabel('x_n')
subplot(2,4,[5 6]);
plot(y); title(['time series y, var:' num2str(var(y))])
xlabel('n'); ylabel('y_n')

subplot(2,4,3);
plot(f,Px,'k')
xlim([0 fs/2]); title(['PSD x, P_{TOT}=' num2str(Ptotx)])
xlabel('f'); ylabel('P_X(f)')
subplot(2,4,7);
plot(f,Py)
xlim([0 fs/2]); title(['PSD y, P_{TOT}=' num2str(Ptoty)])
xlabel('f'); ylabel('P_Y(f)')

subplot(2,4,[4 8]);
plot(ft,Cxyt,'k--','LineWidth',1.2);
hold on; 
plot(f,Cxy,'r','LineWidth',1.2)
xlim([0 fs/2]); ylim([0 1]); title('Coh x,y')
xlabel('f'); ylabel('C^2_{XY}(f)')
legend('true Coh','estimated Coh')



