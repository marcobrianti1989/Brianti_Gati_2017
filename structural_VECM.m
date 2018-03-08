function [B, Xi] = structural_VECM(alp,bet,Gam,res,sigma,nlags,r)

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
T = size(res,1);

%Building matrix Xi - Useful for the long run relations
G = eye(nvar); % Need to build G with a loop because 
%it is (nvar x nvar) while beta is ((nlags*nvar), nvar)
Gam = Gam'; % I do this because I will fill up G always
%subtractin blocks of B corresponding to each lag
for i = 1:nlags
      G = G - Gam(:,(i-1)*nvar+1:i*nvar);
end
Xi = bet*(alp'*G*bet)^(-1)*alp';


obj = @(B) LL_VECM(T,B,sigma); %sum(sum(T/2*logm(B.^2) + T/2*trace((B')^(1)*B^(1)*sigma)));

%Constraint that r(r-1)/2 elements of B must be zero.
Me       = 1; %  Me = no. of equality constraints
Beq      = zeros(Me,1); % Beq is (Me x 1) where
Aeq      = zeros(Me,nvar*nvar); % Aeq is (Me x (nvar*nvar)) - it is B vectorized
Aeq(1,1) = 1; %zero-impact of news on TFP

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Minimization 
B_zero = rand(nvar);
B_opt = fmincon(obj, B_zero,[],[],Aeq,Beq,[],[],@(B) constraint_long(B,Xi,r,nvar),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)






end