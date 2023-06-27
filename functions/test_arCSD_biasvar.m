%% compute empirical confidence intervals of WC estimator of the PSD of a simulated process
clear; close all; clc

%%% Parameters
N=300; % length of simulated time series (increase to see well autocorrelation)
numsimu=200;%quante volte riproduco il processo per stimare bias e varianza

%%%% parameters
SelCrit='BIC'; % model order selection criterion: AIC or BIC
pmax=10; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)
fs=1;

%% generate simu and compute PSD
Ms=2;
c=1.0; %coupling from respiration to RR
par.poles{1}=([0.9 0.25]); % Oscillations resp
par.poles{2}=([0.8 0.1]); % Oscillations RR
par.coup=[1 2 2 c]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[2 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(Ms,par); % parameters
for n=1:numsimu
    clc; disp(['simulation ' int2str(n) ' of ' int2str(numsimu)])
    
    % Realization of the simulation
    U = mvnrnd(zeros(1,Ms),Su,N);
    X=generateVAR(Am,U);
    
    %%% Estimated CSD and Coh
    ret=VARorder(X,pmax); % model order selection
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
    MSC=abs(out_e.COH).^2; %coherence
    eCxy(:,n)=squeeze(MSC(1,2,:));
end

eCxym=mean(eCxy')'; % mean PSD
eCxysd=std(eCxy')'; % std. dev. PSD

% spectral analysis and true spectrum
out_t = PSDvar(Am,Su,ceil(nfft),fs);
MSCt=abs(out_t.COH).^2; %coherence
Cxy=squeeze(MSCt(1,2,:)); f=out_t.f;



%% plots and disps
colm=[86 86 158]/255; % color of mean value
colsh=[189 185 219]/255; %color of shades

figure(1);
hold on;
plot(f,Cxy,'k--','LineWidth',1.5);
plot(f,eCxym,'Color',colm,'LineWidth',1.5);
h1=area(f,eCxym+eCxysd);
h2=area(f,eCxym-eCxysd);
plot(f,Cxy,'k--','LineWidth',1.5);
plot(f,eCxym,'Color',colm,'LineWidth',1.5);
legend('true Coherence','est.Coherence, mean','est.PSD, +-SD')
set(h2,'FaceColor','w');
set(h1,'FaceColor',colsh);
set(h1,'EdgeColor',colsh); set(h2,'EdgeColor','w');
xlim([0 fs/2])
xlabel('f')
ylabel('C^2_{xy}(f)')
