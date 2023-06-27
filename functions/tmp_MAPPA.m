%% calculation of the PSD of a representative EEG
clear; close all; clc;

subjfold='01'; % subject (name of folder) - from 01 to 21
condition='aperti'; % 'aperti' or 'chiusi'

electrode=19;  % EEG to show in first figure: from 1 to 19


%% load EEG data and show one signal
directory = '..\data\';
load([directory subjfold '\' condition '\sottocamp_data_viever.mat']);
data=v_salvato.ser;
chan=v_salvato.nomec';
fs=v_salvato.fc;

data_ch = data(:,11+electrode);
ch_name = chan(11+electrode);

% figure: show EEG of a pre-selected electrode
N=length(data_ch);
t=[1/fs:1/fs:N/fs]'; % time axis
figure(1)
titolo=['Electrode ', char(ch_name)];
title(titolo);
plot(t,data_ch,'b');
xlim([t(1) t(N)])
xlabel('time [sec]');
text(t(10),0.95*max(data_ch),['variance = ' num2str(var(data_ch))])
title(titolo)


%% MAP AN INDEX (variance of EEG time series) across the scalp
%%%% create map (define X and Y coordinates for each electrode)
posX(1) = 3.8;posY(1) = 9.5;
posX(2) = 6.7;posY(2) = 9.5;
posX(3) = 1.2;posY(3) = 7.8;
posX(4) = 3.5;posY(4) = 7.5;
posX(5) = 5.0;posY(5) = 7.5;
posX(6) = 7.0;posY(6) = 7.5;
posX(7) = 8.8;posY(7) = 7.8;
posX(8) = 0.5;posY(8) = 5.0;
posX(9) = 2.8;posY(9) = 5.0;
posX(10) = 5.0;posY(10) = 5.0;
posX(11) = 7.4;posY(11) = 5.0;
posX(12) = 9.5;posY(12) = 5.0;
posX(13) = 1.2;posY(13) = 2.1;
posX(14) = 3.5;posY(14) = 2.5;
posX(15) = 5.0;posY(15) = 2.5;
posX(16) = 7.0;posY(16) = 2.5;
posX(17) = 8.8;posY(17) = 2.1;
posX(18) = 3.8;posY(18) = 0.5;
posX(19) = 6.7;posY(19) = 0.5;

posXn =posX./10; %normalized positions
posYn =posY./10;

K=100; %dimension of the image (in pixels)
[Xi,Yi]=meshgrid(1:1:K,1:K);% create grid (X and Y coordinates of each pixel

%%%% compute index (EEG variance) and interpolate its values on map positions
num_el=19;
for i_el=1:num_el
    data_ch = data(:,11+i_el);
    indice(i_el,:)=var(data_ch);
end
mappa=griddata(posXn*K,posYn*K,indice,Xi,Yi,'v4'); % interpolation of values of indice

%%%% figure: plots interpolated index and shows position and name of electrodes
figure(2); 
hold on; colormap(gca,'jet')
imagesc([0 K],[0 K],mappa);
xlim([0 K]); ylim([0 K]);
colorbar
caxis([min(indice) max(indice)]); %colormap limits
pbaspect([1 1 1]); %square image
set(gca,'FontSize',12,'YTickLabel',[],'XTickLabel',[]); %removes axis labels
% add electrodes and names
for i_el=1:num_el
    plot(posXn(i_el)*K,posYn(i_el)*K,'ok','MarkerFaceColor', [0.5 0.5 0.5],'MarkerSize',18);
    text(posXn(i_el)*K,posYn(i_el)*K,chan(11+i_el),'HorizontalAlignment','center')
end
title('EEG variance');

%%% add circle to the map (simulate scalp)
th = 0:pi/50:2*pi; x= K/2; y= K/2;
xunit = K/2 * cos(th) + x;
yunit = K/2 * sin(th) + y;
plot(xunit, yunit,':k','LineWidth',3);
