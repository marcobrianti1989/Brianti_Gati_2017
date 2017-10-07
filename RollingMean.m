function Y=RollingMean(X,K)
%
% Luca Benati
% European Central Bank
% Monetary Policy Strategy Division
% March 2010
%
[T,N]=size(X);
Y=[];
for tt=K:T
    Y=[Y; mean(X(tt-K+1:tt,:))];
end


















