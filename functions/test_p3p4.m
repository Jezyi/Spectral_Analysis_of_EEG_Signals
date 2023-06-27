% p4 n.27 e p3 n.25
%occhi chiusi Palpha

%% calculation of the PSD of the EEG for all subjects and chamnnels
clear; close all; clc;

condition{1}='aperti'; condition{2}='chiusi'; 

falpha=[8 12]; % limits of ahlpha band (Hz)

ic = 2 %occhi chiusi

%%%% PSD parameters
winname='blackman'; 
typecorrest='unbiased'; 
M=120; % lag at which to truncate the correlation
nfft=1000; %number of points on frequency axis (total)

%% load EEG data + spectral analysis
directory = '..\data\';

    for is=1:21 % 21 subfolders, one for each subject
        clc; disp(['subject ' int2str(is) ', ' condition{ic}])

        % determine string corresponding to subfolder
        if is<10
            subjfold=strcat('0',int2str(is));
        else
            subjfold=int2str(is);
        end

        load([directory subjfold '\' condition{ic} '\sottocamp_data_viever.mat']);
        data=v_salvato.ser;
        chan=v_salvato.nomec';
        fs=v_salvato.fc;

        
        el = [25 27]; %posizione di P3 e P4 
        for i_el=1:length(el)
            x = data(:,el(i_el));

            [Px(:,i_el),f,Pxtot(:,i_el),ftot,Omega]=PSDwc(x,M,typecorrest,winname,fs,nfft);
            
        end

        end
%% statistics
%%% paired t-tests
[h1,p1]=ttest(Px(:,1),Px(:,2)); % aperti vs. chiusi
disp(['h1 = ' num2str(h1)])
disp(['p1 = ' num2str(p1)])
