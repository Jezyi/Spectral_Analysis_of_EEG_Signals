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
        for 
        x = data(:,11+electrode);
        
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



%% plot mean and std of the PSD across all subjects, in the 2 conditions, for one selected channel
Pxm1=mean(Pxv{1},2); % mean PSD, aperti
Pxsd1=std(Pxv{1},0,2)/sqrt(21); % std.error PSD, aperti
Pxm2=mean(Pxv{2},2); % mean PSD, chiusi
Pxsd2=std(Pxv{2},0,2)/sqrt(21); % std.error PSD, chiusi

colm1=[86 86 158]/255; % color of mean value
colsh1=[189 185 219]/255; %color of shades
colm2=[150 0 0]/255; % color of mean value
colsh2=[255 150 150]/255; %color of shades

figure(1);
subplot(2,1,1)
hold on;
ht1=area(f,Pxm1+Pxsd1);
ht2=area(f,Pxm1-Pxsd1);
plot(f,Pxm1,'Color',colm1,'LineWidth',1.5);
set(ht2,'FaceColor','w');
set(ht1,'FaceColor',colsh1);
set(ht1,'EdgeColor',colsh1); set(ht2,'EdgeColor','w');
xlim([0 30])
ylim([0 max(Pxm1+Pxsd1)])
xlabel('f')
ylabel('P_x(f)')
title(['PSD (mean+-std.err), electrode ' chan{11+electrode} ', cond: ' condition{1}])

subplot(2,1,2)
hold on;
ht1=area(f,Pxm2+Pxsd2);
ht2=area(f,Pxm2-Pxsd2);
plot(f,Pxm2,'Color',colm2,'LineWidth',1.5);
set(ht2,'FaceColor','w');
set(ht1,'FaceColor',colsh2);
set(ht1,'EdgeColor',colsh2); set(ht2,'EdgeColor','w');
xlim([0 30])
ylim([0 max(Pxm2+Pxsd2)])
xlabel('f')
ylabel('P_x(f)')
title(['PSD (mean+-std.err), electrode ' chan{11+electrode} ', cond: ' condition{2}])

%% statistics
%%% paired t-tests
[h1,p1]=ttest(Ptot(:,1),Ptot(:,2)); % aperti vs. chiusi
[h2,p2]=ttest(Palpha(:,1),Palpha(:,2)); % aperti vs. chiusi
[h3,p3]=ttest(Palpha_n(:,1),Palpha_n(:,2)); % aperti vs. chiusi

%%% boxplots of the power measures for the two conditions
figure(2)
subplot(1,3,1);
ymin=min([min(Ptot) min(Ptot)]);
ymax=max([max(Ptot) max(Ptot)]);
dy=(ymax-ymin)/20;
boxplot(Ptot,'Widths',0.3,'colors',[0 0.5 1],'symbol','.','positions',[1 2])
ylim([ymin-dy ymax+2*dy])
title('total EEG power')
xticklabels(condition);
ylabel('P_{TOT}')

subplot(1,3,2);
ymin=min([min(Palpha) min(Palpha)]);
ymax=max([max(Palpha) max(Palpha)]);
dy=(ymax-ymin)/20;
boxplot(Palpha,'Widths',0.3,'colors',[0.2 0.2 0.2],'symbol','.','positions',[1 2])
ylim([ymin-dy ymax+2*dy])
title('alpha power')
xticklabels(condition);
ylabel('P_{\alpha}')

subplot(1,3,3);
ymin=min([min(Palpha_n) min(Palpha_n)]);
ymax=max([max(Palpha_n) max(Palpha_n)]);
dy=(ymax-ymin)/20;
boxplot(Palpha_n,'Widths',0.3,'colors',[0.5 0.5 0.5],'symbol','.','positions',[1 2])
ylim([ymin-dy ymax+2*dy])
title('alpha power, %')
xticklabels(condition);
ylabel('P_{\alpha}(%)')


%%% display results of statistical analysis
clc;
disp(['analyzed electrode: ' chan{11+electrode}])
if h1==1, disp(['aperti vs. chiusi: Ptot=' num2str(p1)]); end
if h2==1, disp(['aperti vs. chiusi: Palpha=' num2str(p2)]); end
if h3==1, disp(['aperti vs. chiusi: Palpha%=' num2str(p3)]); end





