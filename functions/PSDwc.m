%% Power Spectral Density with Weighted covariance estimator

function [Pxpos,fpos,Px,f,Omega,rx,rxw,wi,lags,nfft]=PSDwc(x,M,typecorrest,winname,fs,nfft)

narginchk(1,6);% from 1 to 6 input arguments
if nargin < 6, nfft=1001; end 
if nargin < 5, fs=1; end 
if nargin < 4, winname='bartlett'; end 
if nargin < 3, typecorrest='biased'; end 
if nargin < 2, M=length(x); end 

if nfft/2==fix(nfft/2), nfft=nfft+1; end %if nfft is even, make it odd (to compute the PSD at zero lag)

x=x-mean(x); %always work with zero-mean data

[rx,lags]=xcorr(x,M-1,typecorrest); % correlation estimate, truncated at lag M

eval(strcat('wi=',winname,'(',int2str(2*M-1),');'));
rxw=rx.*wi; % truncated and windowed correlation

%%%% DTFT
Omega=(-pi:2*pi/(nfft-1):pi)';
% % Px=nan*ones(nfft,1);
% % for io=1:nfft
% %     Px(io)=exp(-1i*Omega(io)*lags)*rxw;
% % end
Px=exp(-1i*Omega*lags)*rxw;
Px=abs(Px);

f=(-fs/2:fs/((nfft-1)):fs/2)'; % frequency axis

fpos=f(floor(nfft/2)+1:nfft);
Pxpos=Px(floor(nfft/2)+1:nfft);

% % Px2=abs(fft(rxw,nfft)); % verification: fft estimate
% % Px2pos=Px2(floor(nfft/2)+1:nfft);
% % Px=(1/length(x))*(abs(fft(x,nfft)).^2); % verification: periodogram




end
