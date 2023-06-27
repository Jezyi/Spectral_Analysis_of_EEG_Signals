%% test of CSD estimators on simulated RR-RESP time series
clear; close all; clc;

fs=1; % sampling frequency corresponding to mean RR=1 sec
N=300; % length of simulated time series

%%%% parameters
SelCrit='BIC'; % model order selection criterion: AIC or BIC
pmax=10; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)

%% Simulation of cardiorespiratory interactions
% p=5; % model order
Ms=2;
c=0.2; %coupling from respiration to RR
par.poles{1}=([0.9 0.25]); % Oscillations resp
par.poles{2}=([0.8 0.1]); % Oscillations RR
par.coup=[1 2 2 c]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[2 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(Ms,par); % parameters

% Realization of the simulation
U = mvnrnd(zeros(1,Ms),Su,N);
X=generateVAR(Am,U);

%% Cross Spectral Density and Coherence
% cross-spectral analysis and true spectrum
out_t = PSDvar(Am,Su,nfft,fs);
MSCt=abs(out_t.COH).^2; %coherence
Cxyt=squeeze(MSCt(1,2,:)); f=out_t.f;

%%% Estimated CSD and Coh
% model order selection
ret=VARorder(X,pmax);
switch SelCrit
    case 'AIC'
        ep=ret.pottaic;
    case 'BIC'
        ep=ret.pottbic;
end

% VAR model identification
jv=[1 2]; % index of predicted series
iv=[1 2]; % indexes of predictors
iv_lags=(1:ep); % lags of predictors
out_e=LinReg(X,jv,iv,iv_lags);
eAm=out_e.eA;
eSu=out_e.es2u;

% spectrum estimation
out_e = PSDvar(eAm,eSu,nfft,fs);
Px=squeeze(abs(out_e.P(1,1,:)));
Py=squeeze(abs(out_e.P(2,2,:)));
MSC=abs(out_e.COH).^2; %coherence
Cxy=squeeze(MSC(1,2,:));

%%%%%%%%%%%%%%% verifiche con miei codici %%%%%%%%%%%%%%%
% [eAm2,eSu2]=idMVAR(X',ep,0);
% % [[Am; Su] [eAm; eSu] [eAm2'; eSu]]
% [DC2,DTF2,PDC2,GPDC2,COH2,PCOH2,PCOH22,H2,S2,P2,f2]=fdMVAR(eAm2,eSu2,nfft,fs);
% figure(7); plot(f,Cxyt); hold on; plot(f,Cxy,'k'); plot(f,squeeze(abs(COH2(1,2,:)).^2),'r--'); 
%%%%%%%%%%%%%%% verifiche con miei codici %%%%%%%%%%%%%%%

% to get the integral, each PSD point is multiplied for 2*pi/nfft; the integral is then divided by 2*pi to get the variance (see slide 20)
Pxtot=sum(Px)/nfft; Pytot=sum(Py)/nfft;

%% plots
x=X(:,1); y=X(:,2);

figure(1);
subplot(2,4,[1 2]);
plot(x,'k'); title(['time series x, var:' num2str(var(x))])
xlabel('n'); ylabel('x_n')
subplot(2,4,[5 6]);
plot(y); title(['time series y, var:' num2str(var(y))])
xlabel('n'); ylabel('y_n')

subplot(2,4,3);
plot(f,Px,'k')
xlim([0 fs/2]); title(['PSD x, P_{TOT}=' num2str(Pxtot)])
xlabel('f'); ylabel('P_X(f)')
subplot(2,4,7);
plot(f,Py)
xlim([0 fs/2]); title(['PSD y, P_{TOT}=' num2str(Pytot)])
xlabel('f'); ylabel('P_Y(f)')

subplot(2,4,[4 8]);
plot(f,Cxyt,'k--','LineWidth',1.2);
hold on; 
plot(f,Cxy,'r','LineWidth',1.2)
xlim([0 fs/2]); ylim([0 1]); title('Coh x,y')
xlabel('f'); ylabel('C^2_{XY}(f)')
legend('true Coh','estimated Coh')



