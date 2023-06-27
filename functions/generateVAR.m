%% GENERATES A VARPROCESS BY FILTERING A VECTOR NOISE INPUT SERIES
%% X(n)=X(n-1)A(1)+...+X(n-p)A(p)+U(n)

%%% INPUT
% A=[A(1)...A(p)]: M*pM matrix of the MVAR model coefficients (strictly causal model)
% U: N*M matrix of innovations

%%% OUTPUT
% X: N*M matrix of simulated time series

function X=generateVAR(A,U)

[N,M]=size(U);
p=size(A,1)/M;

% X(n)=X(n-1)A(1)+...+X(n-p)A(p)+U(n)
X=zeros(N,M);
for n=1:N
    for k=1:p
        if n-k<=0, break; end % if n<=p, stop when k>=n
        X(n,:)=X(n,:) + X(n-k,:) * A((k-1)*M+(1:M),:) ;
    end
    X(n,:)=X(n,:)+U(n,:);
end


end
