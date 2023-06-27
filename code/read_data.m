clear; close all; clc;

n = 21;    %n. soggetto

condition='chiusi';  % experimental condition: open, closed
electrode=11;  % which EEG channel to open; from 1 to 19


%% load EEG data
for i=1:n

    if i < 10
        subject = ['0' int2str(n)];
    else
        subject = [int2str(n)];
    end
end
%if n<10, subject=['0' int2str(n)];
%else, subject = int2str(n);
%end

percorso=['data\' subject '\' condition '\'];
nome='sottocamp_data_viever';
estensione='.mat';

load([percorso nome estensione]);
data=v_salvato.ser;
chan=v_salvato.nomec';
fs=v_salvato.fc;

% extract EEG signals from the matrix
Y=data(:,12:30); %EEGs are from channel 12 to channel 30 (19 signals)
Y_sel=Y(:,i_el);
Ysel_label = chan(11+i_el); %EEG channel name


%% re-referencing: Common Average Reference
tmp=Y; tmp(:,i_el)=[];
avgY=mean(tmp,2);
Y_sel_reref=Y_sel-avgY;


%% show EEGs 
[N,M]=size(Y);
t=[1/fs:1/fs:N/fs]'; % time axis

%%% show originally acquired data
figure(1)
subplot(2,1,1)
titolo=['subject: ', subject, ' - channel ', char(Ysel_label), ' - condition: ' ,condition , '- raw data'];
title(titolo);
plot(t,Y_sel,'b','LineWidth',1);
xlim([0 max(t)])
xlabel('time (s)');
set(gca,'fontsize',12)
title(titolo)

%%% show re-referenced data 
subplot(2,1,2)
titolo=['subject: ', subject, ' - channel ', char(Ysel_label), ' - condition: ' ,condition , '- re-ref'];
plot(t,Y_sel_reref,'r','LineWidth',1);
xlim([0 max(t)])
xlabel('time (s)');
title(titolo)
set(gca,'fontsize',12)


%% exe
% cfr. referencing: figure(2); plot(t,Y_sel); hold on;  plot(t,Y_sel_reref,'r');
% map: will be seen in exe

%%% compare eyes open and close in occipital channel:
% percorso1=['data\' subject '\aperti\']; percorso2=['data\' subject '\chiusi\'];
% a1=load([percorso1 nome estensione]); a2=load([percorso2 nome estensione]);
% i=29; X1=a1.v_salvato.ser(:,i); X2=a2.v_salvato.ser(:,i); can_i=can{i};
% figure(2); plot(X1); hold on; plot(X2,'r'); title(can_i); legend('open','closed')
% disp([std(X1) std(X2)])

