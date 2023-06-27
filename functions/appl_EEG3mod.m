%% calculation of the PSD of the EEG for all subjects and chamnnels
clear; close all; clc;

condition{1}='aperti'; condition{2}='chiusi'; 

falpha=[8 12]; % limits of ahlpha band (Hz)

%%%% PSD parameters
winname='blackman'; 
typecorrest='unbiased'; 
M=120; % lag at which to truncate the correlation
nfft=1000; %number of points on frequency axis (total)

%% load EEG data + spectral analysis
directory = '..\data\';
for ic =1:2
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

        num_el=19;
        for i_el=1:num_el
            x = data(:,11+i_el);

            [Px,f,Pxtot,ftot,Omega]=PSDwc(x,M,typecorrest,winname,fs,nfft);
            Nf=length(Px); %number of frequency points
            %Pxv{ic,i_el}(:,is)=Px; % export all PSD
            Ptot{ic}(is,i_el)=sum(Px)/Nf;

            % f:n=(fs/2):Nf
            nalpha=round(falpha*Nf*2/fs); % range alpha in frequency bins
            if nalpha(1)==0, nalpha(1)=1; end %adjust inferior bound if needed
            if nalpha(2)==Nf, nalpha(2)=Nf; end %adjust superior bound if needed
            Palpha{ic}(is,i_el)=sum(Px(nalpha(1):nalpha(2)))/Nf;
            Palpha_n{ic}(is,i_el)=Palpha{ic}(is,i_el)/Ptot{ic}(is,i_el);

        end

    end
    
    % compute mean power values for each electrode
    Ptot_m(:,ic)=mean(Ptot{ic})'; % mean across subjects, transpose
    Palpha_m(:,ic)=mean(Palpha{ic})'; % mean across subjects, transpose
    Palpha_n_m(:,ic)=mean(Palpha_n{ic})'; % mean across subjects, transpose
    
end

%% MAP the power across the scalp
% choose index
istring='Palpha_n_m'; % 'Ptot_m' , 'Palpha_m' , 'Palpha_n_m'
eval(['indice1=' istring '(:,1);']); %aperti
eval(['indice2=' istring '(:,2);']); %aperti

% common minimum and maximum values for colormaps
indicemin=min([indice1; indice2]);
indicemax=max([indice1; indice2]);

figure(1); 
subplot(2,1,1)
K=100;
drawmap(indice1,indicemin,indicemax,chan,K)
title([istring 'eyes open']);

subplot(2,1,2)
drawmap(indice2,indicemin,indicemax,chan,K)
title([istring 'eyes closed']);


