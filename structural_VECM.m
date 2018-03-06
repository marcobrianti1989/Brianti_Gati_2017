function [B, Xi] = structural_VECM(alp,bet,Gam,res,sigma,nlags)

% INPUTS:
% res is (T, nvar) - residuals from the reduced form
% sigma is variance-covariance matrix of res
% Gam is (nvar*nlags,nvar) and is the "beta" for the lagged differentiated terms
% Gam does not allowed for the constant so far!
% bet is the cointegration matrix (nvar,r)
% alp is the loading matrix (nvar,r) where r is the number of temporary shocks

% OUTPUT:
% TBD

%Again I am using the notation of Lutkepohl (2005) - EUI, working paper

%Useful tools
nvar = size(sigma,1);

%Building matrix Xi - Useful for the long run relations
G = eye(nvar); % Need to build G with a loop because 
%it is (nvar x nvar) while beta is ((nlags*nvar), nvar)
Gam = Gam'; % I do this because I will fill up G always
%subtractin blocks of B corresponding to each lag
for i = 1:nlags
      G = G - Gam(:,(i-1)*nvar+1:i*nvar);
end
Xi = bet*(alp'*G*bet)^(-1)*alp';

B = 0;






end