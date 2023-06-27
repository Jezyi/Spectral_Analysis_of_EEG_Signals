%% FREQUENCY DOMAIN VAR ANALYSIS from VAR PARAMETERS

%%% inputs: 
% Am=[A(1);...;A(p)]: pM*M matrix of the VAR model coefficients (strictly causal model)
% Su: M*M covariance matrix of the input noises
% nfft= number of points for calculation of the spectral functions
% Fs= sampling frequency

%%% outputs:
% P= Spectral Matrix
% S= Inverse Spectral Matrix
% DC= Directed Coherence
% PDC= Partial Directed Coherence
% COH= Coherence
% PCOH= Partial Coherence
% H= Tranfer Function Matrix
% f= frequency vector

function out = PSDvar(Am,Su,nfft,Fs)
    
    Am=Am'; % the function works with coefficients in row

    M= size(Am,1); % Am has dim M*pM
    p = size(Am,2)/M; % p is the order of the MVAR model

    if nargin<2, Su = eye(M,M); end % if not specified, we assume uncorrelated noises with unit variance as inputs 
    if nargin<3, nfft = 512; end
    if nargin<4, Fs= 1; end    
    if all(size(nfft)==1)	 %if N is scalar
        f = (0:nfft-1)*(Fs/(2*nfft)); % frequency axis
    else            % if N is a vector, we assume that it is the vector of the frequencies
        f = nfft; nfft = length(nfft);
    end

    % s = exp(1i*2*pi*f/Fs); % vector of complex exponentials
    z = 1i*2*pi/Fs;


    %% Initializations: spectral matrices have M rows, M columns and are calculated at each of the N frequencies
    H=zeros(M,M,nfft); % Transfer Matrix
    P=zeros(M,M,nfft); % Spectral Matrix
    S=zeros(M,M,nfft); % Inverse Spectral Matrix
    COH=zeros(M,M,nfft); %Coherence
    PCOH=zeros(M,M,nfft); %Partial Coherence - defined as Dahlhaus 2000
    DC=zeros(M,M,nfft); % directed coherence - defined as Baccala 1998
    PDC=zeros(M,M,nfft); % PDC Baccala Sameshima 2001
    tmp1=zeros(M,1); %denominator for DC (column!)
    tmp2=tmp1; %denominator for DC (column!)
    tmp4=tmp1'; %denominators for PDC (row!)

    A = [eye(M) -Am]; % matrix from which M*M blocks are selected to calculate spectral functions
%     invSu=inv(Su);

    % I define the following matrices forced to be diagonal even when the original Su is not diagonal (this because DC and PDC/GPDC do not use off-diag terms)
    Cd=diag(diag(Su)); % Cd is useful for calculation of DC
    invCd=inv(Cd);% invCd is useful for calculation of GPDC
    %note: in the definition of the DC here (without inst.eff.) the denominator is not the spectrum because the actual Su is not diagonal

    %% computation of spectral functions
    for n=1:nfft % at each frequency

            %%% Coefficient matrix in the frequency domain
            As = zeros(M,M); % matrix As(z)=I-sum(A(k))
            for k = 1:p+1
                As = As + A(:,k*M+(1-M:0))*exp(-z*(k-1)*f(n));  %indicization (:,k*M+(1-M:0)) extracts the k-th M*M block from the matrix B (A(1) is in the second block, and so on)
            end

            %%% Transfer matrix
            H(:,:,n)  = inv(As);

            %%% Spectral matrix
            P(:,:,n)  = H(:,:,n)*Su*H(:,:,n)'; % ' stands for Hermitian transpose

            %%% Inverse Spectral matrix
            S(:,:,n) = inv(P(:,:,n)); % P1(:,:,n) = As'*invSu*As;

            %%% denominators of DC, PDC, GPDC for each m=1,...,num. channels
            for m = 1:M
                tmp1(m)=sqrt((abs(H(m,:,n)).^2) * diag(Cd)); % for the DC: m-th row of H * variance of W (Cd is diag)
                tmp2(m)=sqrt((abs(H(m,:,n)).^2) * ones(M,1)); % for the DTF - don't use covariance information
                tmpp1 = squeeze(As(:,m)); % this takes the m-th column of As...
                tmp4(m) = sqrt(tmpp1'*invCd*tmpp1); % for the GPDC - uses only diagonal covariance information
            end

            %%% Directed Coherence
            DC(:,:,n) = H(:,:,n)*sqrt(Cd) ./ tmp1(:,ones(M,1)); 
            %nota: tmp1(:,ones(M,1)) crea M colonne tutte uguali a tmp1 - la riga (ossia il den) è la stessa - trova in un colpo solo tutti i denominatori per DC

            %%% (Generalized) Partial Directed Coherence
            PDC(:,:,n) = (sqrt(invCd)*As) ./ tmp4(ones(1,M),:);
            
            %NOTE: invertita w.r.t. fdMVAR, perchè anche H lo è (trattazione del corso: il processo è in riga nel VAR)
            PDC(:,:,n)=PDC(:,:,n)';
            DC(:,:,n)=DC(:,:,n)';
    end

    %%% COHERENCE and PARTIAL COHERENCE
    for m1=1:M
        for m2=1:M
            COH(m1,m2,:) = (P(m1,m2,:)) ./ sqrt(abs(P(m1,m1,:).*P(m2,m2,:)));
            PCOH(m1,m2,:) = (-S(m1,m2,:)) ./ sqrt(abs(S(m1,m1,:).*S(m2,m2,:)));      
        end
    end

    out.P=P;
    out.f=f;
    out.H=H;
    out.COH=COH;
    out.PCOH=PCOH;
    out.DC=DC;
    out.PDC=PDC;
    out.P1=S;
end



