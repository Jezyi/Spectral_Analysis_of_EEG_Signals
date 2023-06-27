clear; close all; clc;

condizione{1} = 'aperti'; condizione{2}='chiusi';
falpha = [8 12];                 %limite della banda di frequenza: da 8 a 12 Hertz (ovvero la banda alpha dei segnali alpha)
%parametri del PSD
finestra='blackman';
typecorrest = 'unbiased';        %tipo di correlazione scelta (BT)
M = 120;
nfft = 1000;                     %punti totali nell'asse delle frequenze

%%caricamento dati EEG 
for i = 1:2        %due sottocartelle di 21 sottocartelle :)
    for j = 1:21   %in quanto consideriamo 21 sottocartelle
        clc; disp (['soggetto ' int2str(j) ',' condizione{i}])   %int2str-->trasforma un numero INTero in STRinga
        %scrittura della stringa che determina la sottocartella richiesta
        if j < 10
            sottocartella = strcat ('0',int2str(j));
        else
            sottocartella = int2str(j);
        end
        load (['data\' sottocartella '\' condizione{i} '\sottocamp_data_viever.mat']);
        data = v_salvato.ser;
        chan = v_salvato.nomec;
        fs = v_salvato.fc;
        
        elettrodi_considerati = 19;
        for contatore = 1:elettrodi_considerati
            x = data(:,11 + contatore);          %inserisco i dati di ogni elettrodo di un soggetto in una singola matrice

            % Estrazione dei segnali EEG dalla matrice
            Y=data(:,12:30);                  
            Y_sel=Y(:,contatore);
            Ysel_label = chan(11+contatore);
            tmp=Y; tmp(:,contatore)=[];
            avgY=mean(tmp,2);
            Y_sel_reref=Y_sel-avgY;
            [N,M]=size(Y);
            t=[1/fs:1/fs:N/fs]'; % time axis
            
            %%% Raccolta dei dati di ogni singolo soggetto nelle sue configurazioni
            figure(j)
            subplot(2,1,1)
            titolo=['soggetto: ', int2str(j), ' - canali ', char(Ysel_label), ' - condizione: ' ,condizione{i} , '- raw data'];
            title(titolo);
            plot(t,Y_sel,'b','LineWidth',1);
            xlim([0 max(t)])
            xlabel('time (s)');
            set(gca,'fontsize',12)
            title(titolo)
            %Analisi spettrale tramite la funzione PSDwc
            [Px,f,Pxtot,ftot,Omega] = PSDwc(x,M,typecorrest,finestra,fs,nfft);
            Nf = length(Px);                           %numeri di punti di frequenza
            Ptot{i}(j,contatore) = sum(Px)/Nf;         %calcolo della potenza totale del PSD
            nalpha = round(falpha*Nf*2/fs);            %range di alpha in frequency bins
            if nalpha(1) == 0,  nalpha(1) = 1;  end    %corregge i limiti inferiori se necessario
            if nalpha(2) == Nf, nalpha(2) = Nf; end     
            Palpha{i}(j,contatore) = sum(Px(nalpha(1):nalpha(2)))/Nf;
            Palpha_n{i}(j,contatore) = Palpha{i}(j,contatore)/Ptot{i}(j,contatore);
        end        
    end %fine dei 21 soggetti in una configurazione
    %Calcolo del valore medio di potenza per ogni elettrodo
    %l'apostrofo indica il trasposto
    Ptot_m(:,i) = mean(Ptot{i})';     %media trasposta tra i soggetti  
    Palpha_m(:,i) = mean(Palpha{i})';
    Palpha_n_m(:,i) = mean(Palpha_n{i})';
end     

% Mappatura di potenza sugli scalpi
istring = 'Palpha_m';   %'Ptot_m' , 'Palpha_m'
eval(['indice1=' istring '(:,1);']); %aperti
eval(['indice2=' istring '(:,2);']); %chiusi


% valori comuni di massimo e minimo per la mappa a colori
indicemin = min([indice1; indice2]);
indicemax = max([indice1; indice2]);

figure(1);
subplot(2,1,1)
K = 100;
drawmap (indice1,indicemin,indicemax,chan,K)
title([istring 'eyes open']);

subplot (2,1,2)
drawmap (indice2,indicemin,indicemax,chan,K)
title([istring 'eyes closed']);
