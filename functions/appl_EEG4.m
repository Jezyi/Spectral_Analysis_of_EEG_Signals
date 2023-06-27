%% VAR spectral decomposition for three channels - single file analysis
clear; close all; clc;

subjfold='04'; % subject (name of folder) - from 01 to 21
condition='aperti'; % 'aperti' or 'chiusi'
electrodes=[1 9 18];  % front-back Fp1-C3-O1: [1 9 18] , Fp2-C4-O2: [2 11 19]

fband=[8 12]; bandname='alpha'; % limits of selected band (alpha:8-13; beta:13-30)

%%%% Spectral Analysis Parameters
p=6; % fix mode order (selection criteria don't find a minimum for this dataset)
pmax=20; %maximum scanned model order
nfft=1000; %number of points on frequency axis (positive frequencies)


%% load EEG data and show the three signals
directory = '..\exe_01_Database\D-EEG\data\';
load([directory subjfold '\' condition '\sottocamp_data_viever.mat']);
data=v_salvato.ser;
chan=v_salvato.nomec';
fs=v_salvato.fc;

Xo = data(:,11+electrodes);
ch_name = chan(11+electrodes);

[N,M]=size(Xo);
%%% filter to remove slow trends
pfilter=0.95;
for m=1:M % remove mean from all time series
    Xf(:,m)=AR_filter(Xo,m,pfilter); % filtered series
    X(:,m)=Xf(:,m)-mean(Xf(:,m));
end


%% spectral analysis
%%% model order
ret=VARorder(X,pmax);
figure(1);
plot(ret.aic,'b.-')
hold on; plot(ret.bic,'k.-')
legend('AIC','BIC')


% VAR model identification
jv=[1:M]; % index of predicted series
iv=[1:M]; % indexes of predictors
iv_lags=(1:p); % lags of predictors
outR=LinReg(X,jv,iv,iv_lags);
eAm=outR.eA;
eSu=outR.es2u;
% estimation of all frequency domain functions
out = PSDvar(eAm,eSu,nfft,fs);
eP=abs(out.P); % spectral matrix
eCOH=abs(out.COH).^2; % (magnitude squared) coherence
eDC=abs(out.DC).^2; % (magnitude squared) directed coherence
eS=abs(out.P1); % inverse spectral matrix
ePCOH=abs(out.PCOH).^2; % (magnitude squared) partial coherence
ePDC=abs(out.PDC).^2; % (magnitude squared) partial directed coherence
f=out.f; Nf=length(f);




%% plot time series and spectra
t=(1/fs:1/fs:N/fs)'; % time axis
%colors for shades
colsh_gray=[150 150 150]/255; 
colsh_red=[255 150 150]/255; 
colsh_orange=[255 200 120]/255; 
colsh_azure=[150 210 240]/255; 
colsh_green=[100 190 190]/255; 
% colsh1=[189 185 219]/255; 
% colsh2=[255 180 180]/255; %color of shades

% f:n=(fs/2):Nf
nband=round(fband*Nf*2/fs); % range alpha in frequency bins
if nband(1)==0, nband(1)=1; end %adjust inferior bound if needed
if nband(2)==Nf, nband(2)=Nf; end %adjust superior bound if needed

h0=figure('numbertitle','off','name','EEG signals');
h1=figure('numbertitle','off','name','Spectra and Coherence');
h2=figure('numbertitle','off','name','Inverse Spectra and Partial Coherence');
h3=figure('numbertitle','off','name','Directed Coherence');
h4=figure('numbertitle','off','name','Partial Directed Coherence');
eCOH_band=nan*ones(M,M);
ePCOH_band=nan*ones(M,M);
eDC_band=nan*ones(M,M);
ePDC_band=nan*ones(M,M);
for i=1:M
    figure(h0); % time series 
    subplot(M,1,i)
    plot(t,X(:,i),'color',[0 0.4 0.7]);
    xlim([t(1) t(N)])
    xlabel('time [sec]');
    title(['EEG signal Electrode ', char(ch_name{i}) ', ' condition])
    
    for j=1:M
        figure(h1); % spectra and coherence       
        subplot(M,M,(i-1)*M+j);
        if i==j %diagonal plots: PSD
            Ptmp=squeeze(eP(i,j,:)); % PSD to show     
            plot(f,Ptmp,'k','Linewidth',1.1);
            xlim([0 fs/2]); title(['P_{' ch_name{i} '}(f)']);
            hold on; ha1=area(f(nband(1):nband(2)),Ptmp(nband(1):nband(2))); set(ha1,'FaceColor',colsh_gray);           
        else % off-diagonal plots: COHERENCE
            COHtmp=squeeze(eCOH(i,j,:));
            plot(f, COHtmp,'r','Linewidth',1.1);
            xlim([0 fs/2]); ylim([0 1]); title(['|\Gamma_{' ch_name{i} '-' ch_name{j} '}(f)|^2']);
            hold on; ha2=area(f(nband(1):nband(2)),COHtmp(nband(1):nband(2))); set(ha2,'FaceColor',colsh_red);
            eCOH_band(i,j)=mean(COHtmp(nband(1):nband(2)));
            text(f(nband(2)),0.9,['\Gamma_{' int2str(i) int2str(j) '}(\' bandname ')=' num2str(eCOH_band(i,j))])
        end
        
        figure(h2); % inverse spectra and partial coherence       
        subplot(M,M,(i-1)*M+j);
        if i==j %diagonal plots: inverse PSD
            plot(f, squeeze(eS(i,j,:)),'Color',[0.75 0.75 0.75],'Linewidth',1.1);
            xlim([0 fs/2]); title(['S_{' ch_name{i} '}(f)']);
        else % off-diagonal plots: PARTIAL COHERENCE
            PCOHtmp=squeeze(ePCOH(i,j,:));
            plot(f,PCOHtmp,'Color',[1 0.5 0],'Linewidth',1.1);
            xlim([0 fs/2]); ylim([0 1]); title(['|\Pi_{' ch_name{i} '-' ch_name{j} '}(f)|^2']);
            hold on; ha2=area(f(nband(1):nband(2)),PCOHtmp(nband(1):nband(2))); set(ha2,'FaceColor',colsh_orange);
            ePCOH_band(i,j)=mean(PCOHtmp(nband(1):nband(2)));
            text(f(nband(2)),0.9,['\Pi_{' int2str(i) int2str(j) '}(\' bandname ')=' num2str(ePCOH_band(i,j))])
        end
        
        figure(h3); % directed coherence       
        subplot(M,M,(i-1)*M+j);
        DCtmp=squeeze(eDC(i,j,:));
        plot(f,DCtmp,'-','Color',[0 0.5 1],'Linewidth',1.1);
        xlim([0 fs/2]); ylim([0 1]); title(['|\gamma_{' ch_name{i} '\rightarrow' ch_name{j} '}(f)|^2']);
        hold on; ha2=area(f(nband(1):nband(2)),DCtmp(nband(1):nband(2))); set(ha2,'FaceColor',colsh_azure);
        eDC_band(i,j)=mean(DCtmp(nband(1):nband(2)));
        text(f(nband(2)),0.9,['\gamma_{' int2str(i) int2str(j) '}(\' bandname ')=' num2str(eDC_band(i,j))])
        
        figure(h4); % partial directed coherence       
        subplot(M,M,(i-1)*M+j);
        PDCtmp=squeeze(ePDC(i,j,:));
        plot(f,PDCtmp,'-','Color',[0 0.5 0.5],'Linewidth',1.1);
        xlim([0 fs/2]); ylim([0 1]); title(['|\pi_{' ch_name{i} '\rightarrow' ch_name{j} '}(f)|^2']);      
        hold on; ha2=area(f(nband(1):nband(2)),PDCtmp(nband(1):nband(2))); set(ha2,'FaceColor',colsh_green);
        ePDC_band(i,j)=mean(PDCtmp(nband(1):nband(2)));
        text(f(nband(2)),0.9,['\pi_{' int2str(i) int2str(j) '}(\' bandname ')=' num2str(ePDC_band(i,j))])
    end
end

