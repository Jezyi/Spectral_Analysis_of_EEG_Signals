%% show poles and spectra of simulated AR process
clear; close all; clc;

fs=1; % sampling frequency
nfft=1000; %number of points on frequency axis (positive frequencies)

simutype='resp'; % 'resp' or 'RR'

%% Simulation of AR process of order 5 (mimicking RR interval)
M=1;
switch simutype
    case 'resp'
        par.poles{1}=([0.8 0.1]); % Oscillations 1 (resp)
        par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
        par.Su=1; %variance of innovation processes
        name='Simulated respiration';
    case 'RR'
        par.poles{1}=([0.65 0; 0.8 0.1; 0.92 0.25]); % Oscillations RR
%         par.poles{1}=([0 0; 0.7 0.1; 0.92 0.35]); % Oscillations RR
        par.coup=[]; %in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 
        par.Su=1; %variance of innovation processes
        name='Simulated RR interval';
end

%%% Theoretical VAR process
[Am,Su,Ak,z]=theoreticalVAR(1,par); % parameters

poles=z{1};
% cartesian coordinates
a=real(poles);
b=imag(poles);
% polar coordinates
rho=abs(poles);
phi=angle(poles);


%% Theoretical (true) PSD
out = PSDvar(Am,Su,nfft,fs);
Px=squeeze(out.P);
fpos=out.f';

% instructions to draw the unit circle
f=[-flip(fpos,1); fpos(2:length(fpos))];
Omega=2*pi*f/fs;
xc=cos(Omega);
yc=sin(Omega);

figure(1)
subplot(2,1,1)
plot(xc, yc,'k','LineWidth',1.5);
line([-1.25 1.25],[0 0],'color','k');
line([0 0],[-1.25 1.25],'color','k');
hold on; plot(a,b,'o','color',[0 0.2 1],'MarkerSize',6,'MarkerFaceColor',[0 0.5 0.8])
pbaspect([1 1 1]); %square image
text(1.34,0,'Re','HorizontalAlignment','left')
text(0,1.34,'Im','HorizontalAlignment','center')
text(1.03,0.1,'f=0','HorizontalAlignment','left','FontSize',10)
text(-1.03,0.1,'f=fs/2','HorizontalAlignment','right','FontSize',10)
xlim([-1.25 1.25]);ylim([-1.25 1.25])
axis off

subplot(2,1,2)
plot(fpos,Px,'LineWidth',1.5)
xlabel('f'); ylabel('P_X(f)')


