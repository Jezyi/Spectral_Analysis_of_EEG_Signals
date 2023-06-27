%% calculation of the PSD of the EEG for all subjects in a selected channel
clear; close all; clc;

condition{1}='aperti'; condition{2}='chiusi'; 
numsubj=21; % number of subjects in the database
electrode= [25 27];  % EEG to analyze P3 P4
falpha=[8 12]; % limits of ahlpha band (Hz)

%%%% PSD parameters
winname='blackman'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
typecorrest='unbiased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
M=120; % lag at which to truncate the correlation
nfft=1000; %number of points on frequency axis (total)

%% load EEG data + spectral analysis
directory = '..\exe_01_Database\D-EEG\data\';

Palpha=nan*ones(numsubj,2); Palpha_n=Palpha; Ptot=Palpha_n;
for ic =1:2
    for is=1:numsubj % 21 subfolders, one for each subject
        clc; disp(['subject ' int2str(is) ', ' condition{ic}])

        % determine string corresponding to subfolder
        if is<10
            subjfold=strcat('0',int2str(is));
        else
            subjfold=int2str(is);
        end
        
        %%% load data
        load([directory subjfold '\' condition{ic} '\sottocamp_data_viever.mat']);
        data=v_salvato.ser;
        chan=v_salvato.nomec';
        fs=v_salvato.fc;
        for i_el=1:length(el)
            x = data(:,electrode(i_el));

      
        
        %%% spectral analysis
        %indice{ic}(is)=var(x);
        [Px,f,Pxtot,ftot,Omega]=PSDwc(x,M,typecorrest,winname,fs,nfft);
        Nf=length(Px); %number of frequency points
        Pxv{ic}(:,is)=Px; % export all PSD
        Ptot(is,ic)=sum(Px)/Nf;

        % f:n=(fs/2):Nf
        nalpha=round(falpha*Nf*2/fs); % range alpha in frequency bins
        if nalpha(1)==0, nalpha(1)=1; end %adjust inferior bound if needed
        if nalpha(2)==Nf, nalpha(2)=Nf; end %adjust superior bound if needed
        Palpha(is,ic)=sum(Px(nalpha(1):nalpha(2)))/Nf;
        Palpha_n(is,ic)=100*Palpha(is,ic)/Ptot(is,ic);

    end
    end
end





%%% display results of statistical analysis
clc;
disp(['analyzed electrode: ' chan{11+electrode}])
if h1==1, disp(['aperti vs. chiusi: Ptot=' num2str(p1)]); end
if h2==1, disp(['aperti vs. chiusi: Palpha=' num2str(p2)]); end
if h3==1, disp(['aperti vs. chiusi: Palpha%=' num2str(p3)]); end





