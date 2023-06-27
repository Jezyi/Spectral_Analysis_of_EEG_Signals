%% Power Spectral Density with Weighted covariance estimator

function [Cxy,Px1,Py1,Pxy1,f]=COHwc(xy,M,typecorrest,winname,fs,nfft)

narginchk(1,6);% from 1 to 6 input arguments
if nargin < 6, nfft=1001; end 
if nargin < 5, fs=1; end 
if nargin < 4, winname='rectwin'; end 
if nargin < 3, typecorrest='biased'; end 
if nargin < 2, M=length(xy); end 

if nfft/2==fix(nfft/2), nfft=nfft+1; end %if nfft is even, make it odd (to compute the PSD at zero lag)

% clear; close all; clc;
% winname='hann'; % 'rectwin' for rectangular, 'bartlett' for triangular, or window name ('blackman','hann','hamming',...)
% typecorrest='biased'; % 'unbiased' for BT estimator, 'biased' for JW estimator
% M=40; % lag at which to truncate the correlation
% nfft=1000; %number of points on frequency axis (total)
% %data
% data=load('CardiovascularExample.txt');
% SAP=data(:,1);
% RR=data(:,3);
% RESP=data(:,4);
% fs=1000/mean(RR);
% xy=[SAP RR];


% two zero-mean time series
x=xy(:,1);y=xy(:,2);
x=x-mean(x); 
y=y-mean(y);

% correlation estimates, truncated at lag M
[rx,lags]=xcorr(x,M-1,typecorrest); 
ry=xcorr(y,M-1,typecorrest); 
rxy=xcorr(x,y,M-1,typecorrest);

% truncated and windowed correlation
eval(strcat('wi=',winname,'(',int2str(2*M-1),');'));
rxw=rx.*wi; 
ryw=ry.*wi;
rxyw=rxy.*wi;



%%%% DTFTs
Omega=(-pi:2*pi/(nfft-1):pi)';
Px=abs(exp(-1i*Omega*lags)*rxw); % PSD of x
Py=abs(exp(-1i*Omega*lags)*ryw); % PSD of y
Pxy=exp(-1i*Omega*lags)*rxyw; % cross-PSD between x and y


ffull=(-fs/2:fs/((nfft-1)):fs/2)'; % frequency axis
% spectral functions between f=0 and f=fs/2
f=ffull(floor(nfft/2)+1:nfft);
Px1=Px(floor(nfft/2)+1:nfft);
Py1=Py(floor(nfft/2)+1:nfft);
Pxy1=Pxy(floor(nfft/2)+1:nfft);

Cxy=(abs(Pxy1).^2)./(Px1.*Py1); %magnitude squared coherence




end
