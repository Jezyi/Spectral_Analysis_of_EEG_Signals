%% test of Partial and Causal Coherence estimators on simulated RR-SAP-RESP time series
clear; close all; clc;

fs=1; % sampling frequency corresponding to mean RR=1 sec
N=300; % length of simulated time series

%%%% parameters
SelCrit='BIC'; % model order selection criterion: AIC or BIC
pmax=10; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)

%% Simulation of cardiovascular and cardiorespiratory interactions
M=3;
c=0; %coupling from respiration to RR
par.poles{1}=([0.9 0.3]); % Oscillations resp
par.poles{2}=([0.8 0.1]); % Oscillations RR
par.poles{3}=([0.8 0.1]); % Oscillations SAP
par.coup=[1 2 1 c; 1 3 1 1-c; 3 2 2 0.5; 2 3 2 0]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
par.Su=[2 1 1]; %variance of innovation processes

%%% Theoretical VAR process
[Am,Su,Ak]=theoreticalVAR(M,par); % parameters

% Realization of the simulation
U = mvnrnd(zeros(1,M),Su,N);
X=generateVAR(Am,U);

%% frequency domain parametric analysis
%%% TRUE VALUES (from known VAR parameters)
out_t = PSDvar(Am,Su,nfft,fs);
f=out_t.f'; % frequency vector
P=abs(out_t.P); % spectral matrix
COH=abs(out_t.COH).^2; % (magnitude squared) coherence
DC=abs(out_t.DC).^2; % (magnitude squared) directed coherence
H=out_t.H; % transfer matrix
S=abs(out_t.P1); % inverse spectral matrix
PCOH=abs(out_t.PCOH).^2; % (magnitude squared) partial coherence
PDC=abs(out_t.PDC).^2; % (magnitude squared) partial directed coherence

%%% ESTIMATED VALUES (from VAR parameters estimated from data)
% model order selection
ret=VARorder(X,pmax);
switch SelCrit
    case 'AIC'
        ep=ret.pottaic;
    case 'BIC'
        ep=ret.pottbic;
end
% VAR model identification
jv=[1 2 3]; % index of predicted series
iv=[1 2 3]; % indexes of predictors
iv_lags=(1:ep); % lags of predictors
out_e=LinReg(X,jv,iv,iv_lags);
eAm=out_e.eA;
eSu=out_e.es2u;
% estimation of all frequency domain functions
out_e = PSDvar(eAm,eSu,nfft,fs);
eP=abs(out_e.P); % spectral matrix
eCOH=abs(out_e.COH).^2; % (magnitude squared) coherence
eDC=abs(out_e.DC).^2; % (magnitude squared) directed coherence
eH=out_e.H; % transfer matrix
eS=abs(out_e.P1); % inverse spectral matrix
ePCOH=abs(out_e.PCOH).^2; % (magnitude squared) partial coherence
ePDC=abs(out_e.PDC).^2; % (magnitude squared) partial directed coherence



%% plots
h1=figure('numbertitle','off','name','Spectra and Coherence');
h2=figure('numbertitle','off','name','Inverse Spectra and Partial Coherence');
h3=figure('numbertitle','off','name','Directed Coherence');
h4=figure('numbertitle','off','name','Partial Directed Coherence');
for i=1:M
    for j=1:M
        figure(h1); % spectra and coherence       
        subplot(M,M,(i-1)*M+j);
        if i==j %diagonal plots: PSD
            plot(f, squeeze(P(i,j,:)),'k--','Linewidth',1.1);
            hold on; plot(f, squeeze(eP(i,j,:)),'k','Linewidth',1.1);
            xlim([0 fs/2]); title(['P_{X' int2str(i) '}(\Omega)']);
            if i==1, legend('PSD,true','PSD,est.'); end
        else % off-diagonal plots: COHERENCE
            plot(f, squeeze(COH(i,j,:)),'r--','Linewidth',1.1);
            hold on; plot(f, squeeze(eCOH(i,j,:)),'r','Linewidth',1.1);
            xlim([0 fs/2]); ylim([0 1]); title(['|\Gamma_{X_' int2str(i) 'X_' int2str(j) '}(\Omega)|^2']);
            if i==1, legend('COH,true','COH,est.'); end
        end
        
        figure(h2); % inverse spectra and partial coherence       
        subplot(M,M,(i-1)*M+j);
        if i==j %diagonal plots: inverse PSD
            plot(f, squeeze(S(i,j,:)),'--','Color',[0.75 0.75 0.75],'Linewidth',1.1);
            hold on; plot(f, squeeze(eS(i,j,:)),'Color',[0.75 0.75 0.75],'Linewidth',1.1);
            xlim([0 fs/2]); title(['S_{X_' int2str(i) '}(\Omega)']);
            if i==1, legend('IPSD,true','IPSD,est.'); end
        else % off-diagonal plots: PARTIAL COHERENCE
            plot(f, squeeze(PCOH(i,j,:)),'--','Color',[1 0.5 0],'Linewidth',1.1);
            hold on; plot(f, squeeze(ePCOH(i,j,:)),'Color',[1 0.5 0],'Linewidth',1.1);
            xlim([0 fs/2]); ylim([0 1]); title(['|\Pi_{X_' int2str(i) 'X_' int2str(j) '}(\Omega)|^2']);
            if i==1, legend('PCOH,true','PCOH,est.'); end
        end
        
        figure(h3); % directed coherence       
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(DC(i,j,:)),'--','Color',[0 0.5 1],'Linewidth',1.1);
        hold on; plot(f, squeeze(eDC(i,j,:)),'-','Color',[0 0.5 1],'Linewidth',1.1);
        xlim([0 fs/2]); ylim([0 1]); title(['|\gamma_{X_' int2str(i) '\rightarrowX_' int2str(j) '}(\Omega)|^2']);
        if i==1, legend('DC,true','DC,est.'); end
        
        figure(h4); % partial directed coherence       
        subplot(M,M,(i-1)*M+j);
        plot(f, squeeze(PDC(i,j,:)),'--','Color',[0 0.5 0.5],'Linewidth',1.1);
        hold on; plot(f, squeeze(ePDC(i,j,:)),'-','Color',[0 0.5 0.5],'Linewidth',1.1);
        xlim([0 fs/2]); ylim([0 1]); title(['|\pi_{X_' int2str(i) '\rightarrowX_' int2str(j) '}(\Omega)|^2']);
        if i==1, legend('PDC,true','PDC,est.'); end
        
        
    end
end
