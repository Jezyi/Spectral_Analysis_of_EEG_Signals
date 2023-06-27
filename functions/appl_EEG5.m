%% calculation of Spectral decomposition for occipital, central and frontal EEG -all channels and subjects
clear; close all; clc;

condition{1}='aperti'; condition{2}='chiusi'; 
numsubj=21; % number of subjects in the database
electrodes=[1 9 18];  % front-back Fp1-C3-O1: [1 9 18] , Fp2-C4-O2: [2 11 19]

fband=[8 12]; bandname='alpha'; % limits of selected band (alpha:8-13; beta:13-30)

%%%% Spectral Analysis Parameters
p=6; % fix mode order (selection criteria don't find a minimum for this dataset)
nfft=1000; %number of points on frequency axis (positive frequencies)

%% load EEG data + spectral analysis
directory = 'C:\Users\Chiara Galfano\Documents\MATLAB\SABS\data\';
for ic =1:2 % ic=1: eyes open; ic=2: eyes closed
    for is=1:numsubj % 21 subfolders, one for each subject
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
                
        Xo = data(:,11+electrodes);
        ch_name = chan(11+electrodes);
        [N,M]=size(Xo);
        pfilter=0.95;
        for m=1:M % remove mean from all time series
            Xf(:,m)=AR_filter(Xo,m,pfilter); % filtered series
            X(:,m)=Xf(:,m)-mean(Xf(:,m));
        end
        
        % VAR model identification
        jv=[1:M]; % index of predicted series
        iv=[1:M]; % indexes of predictors
        iv_lags=(1:p); % lags of predictors
        outR=LinReg(X,jv,iv,iv_lags);
        eAm=outR.eA;
        eSu=outR.es2u;
        % estimation of all frequency domain functions
        out = PSDvar(eAm,eSu,nfft,fs);
        %eP=abs(out.P); % spectral matrix
        eCOH=abs(out.COH).^2; % (magnitude squared) coherence
        eDC=abs(out.DC).^2; % (magnitude squared) directed coherence
        eS=abs(out.P1); % inverse spectral matrix
        ePCOH=abs(out.PCOH).^2; % (magnitude squared) partial coherence
        ePDC=abs(out.PDC).^2; % (magnitude squared) partial directed coherence
        f=out.f; Nf=length(f);
        
        % f:n=(fs/2):Nf
        nband=round(fband*Nf*2/fs); % range alpha in frequency bins
        if nband(1)==0, nband(1)=1; end %adjust inferior bound if needed
        if nband(2)==Nf, nband(2)=Nf; end %adjust superior bound if needed
        for i=1:M
            for j=1:M
                COHtmp=squeeze(eCOH(i,j,:));
                eCOH_band{ic}{i,j}(is)=mean(COHtmp(nband(1):nband(2)));
                PCOHtmp=squeeze(ePCOH(i,j,:));
                ePCOH_band{ic}{i,j}(is)=mean(PCOHtmp(nband(1):nband(2)));
                DCtmp=squeeze(eDC(i,j,:));
                eDC_band{ic}{i,j}(is)=mean(DCtmp(nband(1):nband(2)));
                PDCtmp=squeeze(ePDC(i,j,:));
                ePDC_band{ic}{i,j}(is)=mean(PDCtmp(nband(1):nband(2)));
            end
        end

    end
   
end

%% collect results
%%% COH, PCOH results
cnt=1; id1=[];
for i=1:M
    for j=i+1:M
        % eyes open
        res_COH_eo(:,cnt)=eCOH_band{1}{i,j}';
        res_PCOH_eo(:,cnt)=ePCOH_band{1}{i,j}';
        % eyes closed
        res_COH_ec(:,cnt)=eCOH_band{2}{i,j}';
        res_PCOH_ec(:,cnt)=ePCOH_band{2}{i,j}';
        id1{cnt}=[ch_name{i} ',' ch_name{j}];
        cnt=cnt+1;
    end
end

%%% DC and PDC results
cnt=1; id2=[];
for i=1:M
    for j=i+1:M
        %%% direction i->j
        % eyes open
        res_DC_eo(:,cnt)=eDC_band{1}{i,j}';
        res_PDC_eo(:,cnt)=ePDC_band{1}{i,j}';
        % eyes closed
        res_DC_ec(:,cnt)=eDC_band{2}{i,j}';
        res_PDC_ec(:,cnt)=ePDC_band{2}{i,j}';
        id2{cnt}=[ch_name{i} '->' ch_name{j}];
        cnt=cnt+1;
        
        %%% direction j->i
        % eyes open
        res_DC_eo(:,cnt)=eDC_band{1}{j,i}';
        res_PDC_eo(:,cnt)=ePDC_band{1}{j,i}';
        % eyes closed
        res_DC_ec(:,cnt)=eDC_band{2}{j,i}';
        res_PDC_ec(:,cnt)=ePDC_band{2}{j,i}';
        id2{cnt}=[ch_name{j} '->' ch_name{i}];
        cnt=cnt+1;
    end
end


%% statistical analyses
for m=1:size(res_COH_eo,2)
    [hvalCOH(m),pvalCOH(m)]=ttest(res_COH_eo(:,m),res_COH_ec(:,m)); % aperti vs. chiusi
    [hvalPCOH(m),pvalPCOH(m)]=ttest(res_PCOH_eo(:,m),res_PCOH_ec(:,m)); % aperti vs. chiusi
end

for m=1:size(res_DC_eo,2)
    [hvalDC(m),pvalDC(m)]=ttest(res_DC_eo(:,m),res_DC_ec(:,m)); % aperti vs. chiusi
    [hvalPDC(m),pvalPDC(m)]=ttest(res_PDC_eo(:,m),res_PDC_ec(:,m)); % aperti vs. chiusi
end

%% BOXPLOTS
col_open=[30 80 160]/255;
col_closed=[150 20 20]/255;
xpos=[1 3 5];

figure(1) %% COH
boxplot(res_COH_eo,'Widths',0.3,'colors',col_open,'symbol','.','positions',xpos)
hold on;
boxplot(res_COH_ec,'Widths',0.3,'colors',col_closed,'symbol','.','positions',xpos+0.5)
title(['COH, mean values in ' bandname ' band'])
xticks(xpos+0.25)
xticklabels(id1);
xlim([0 6.5]);ylim([-0.1 1])
text(3,0.97,'eyes open','color',col_open,'HorizontalAlignment','Center')
text(3,0.93,'eyes closed','color',col_closed,'HorizontalAlignment','Center')
for m=1:size(res_COH_eo,2)
    if hvalCOH(m)==1
        text(xpos(m)+0.25,-0.05,'*','Color','r','FontSize',16,'HorizontalAlignment','Center');
    end
end

figure(2) %% PCOH
boxplot(res_PCOH_eo,'Widths',0.3,'colors',col_open,'symbol','.','positions',xpos)
hold on;
boxplot(res_PCOH_ec,'Widths',0.3,'colors',col_closed,'symbol','.','positions',xpos+0.5)
title(['PCOH, mean values in ' bandname ' band'])
xticks(xpos+0.25)
xticklabels(id1);
xlim([0 6.5]);ylim([-0.1 1])
text(3,0.97,'eyes open','color',col_open,'HorizontalAlignment','Center')
text(3,0.93,'eyes closed','color',col_closed,'HorizontalAlignment','Center')
for m=1:size(res_PCOH_eo,2)
    if hvalPCOH(m)==1
        text(xpos(m)+0.25,-0.05,'*','Color','r','FontSize',16,'HorizontalAlignment','Center');
    end
end


    
xpos=[1 3 5 7 9 11];
figure(3) %% DC
boxplot(res_DC_eo,'Widths',0.3,'colors',col_open,'symbol','.','positions',xpos)
hold on;
boxplot(res_DC_ec,'Widths',0.3,'colors',col_closed,'symbol','.','positions',xpos+0.5)
title(['DC, mean values in ' bandname ' band'])
xticks(xpos+0.25)
xticklabels(id2);
xlim([0 12.5]);ylim([-0.1 1])
text(1,0.97,'eyes open','color',col_open,'HorizontalAlignment','Left')
text(1,0.93,'eyes closed','color',col_closed,'HorizontalAlignment','Left')
for m=1:size(res_DC_eo,2)
    if hvalDC(m)==1
        text(xpos(m)+0.25,-0.05,'*','Color','r','FontSize',16,'HorizontalAlignment','Center');
    end
end

figure(4) %% PDC
boxplot(res_PDC_eo,'Widths',0.3,'colors',col_open,'symbol','.','positions',xpos)
hold on;
boxplot(res_PDC_ec,'Widths',0.3,'colors',col_closed,'symbol','.','positions',xpos+0.5)
title(['PDC, mean values in ' bandname ' band'])
xticks(xpos+0.25)
xticklabels(id2);
xlim([0 12.5]);ylim([-0.1 1])
text(1,0.97,'eyes open','color',col_open,'HorizontalAlignment','Left')
text(1,0.93,'eyes closed','color',col_closed,'HorizontalAlignment','Left')
for m=1:size(res_DC_eo,2)
    if hvalPDC(m)==1
        text(xpos(m)+0.25,-0.05,'*','Color','r','FontSize',16,'HorizontalAlignment','Center');
    end
end




