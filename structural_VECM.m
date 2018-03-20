function [B_opt, Xi, A] = structural_VECM(alp,bet,Gam,res,sigma,nlags,r)

% INPUTS:
% res is (nvar,T) - residuals from the reduced form
% sigma is variance-covariance matrix of res
% Gam is (nvar,nvar*nlags) and is the "beta" for the lagged differentiated terms
% Gam does not allowed for the constant so far!
% bet is the cointegration matrix (nvar,r)
% alp is the loading matrix (nvar,r) where r is the number of temporary shocks

% OUTPUT:
% Obtain the impact matrix B (nvar,nvar). Interpretation of B is pretty
% similar to the impact matrist of a SVAR. The procedure is more difficult
% since we need just r(r-1)/2 constraints (impact restrictions) since r*(nvar - r)
% independent restrictions have been already imposed with the reduced rank regression
% from the reduced form VECM. I still slightly miss the intuition. But the
% idea is that number of cointegrations is related to the number of
% transitory shocks since it implies that in the long run, the effect of a
% transitory shocks can be swiped out and it remains just the common trend.
% A represents the structural relation between the past values and the
% present values. A is (nvar*nlags,nvar). Treat A as you used to treat B for a simple SVAR!

% NOTICE
% Again I am using the notation of Lutkepohl (2005) - EUI, working paper

%Useful tools
nvar = size(sigma,1);
T = size(res,2);

%Building matrix Xi - Useful for the long run relations
G = eye(nvar); % Need to build G with a loop because
%it is (nvar x nvar) while beta is ((nlags*nvar), nvar)
for i = 1:nlags
      G = G - Gam(:,((i-1)*nvar)+1:i*nvar); %(nvar,nvar)
end
% Xi = bet*(alp'*G*bet)^(-1)*alp'; %(nvar,nvar)
betT = null(bet','r');
alpT = null(alp','r');
Xi = betT*(alpT'*G*betT)^(-1)*alpT'; %(nvar,nvar) % This line is based on Chiawa et al. 


% Objective function to minimize - See Lutkepohl (2005) - EUI, working paper
obj = @(B) LL_VECM(T,B,sigma);
%sum(sum(T/2*log(B.^2) + T/2*trace((B')^(-1)*B^(-1)*sigma)));

%Constraint that r(r-1)/2 elements of B must be zero.
if r > 1
      Me       = r*(r-1)/2; %  Me = no. of equality constraints. It must be r*(r - 1)/2 in a VECM
      Beq      = zeros(Me,1); % Beq is (Me x 1) where
      Aeq      = zeros(Me,nvar*nvar); % Aeq is (Me x (nvar*nvar)) - it is B vectorized
      element_restricted = 2; %Coordinate to set the element equal to zero. It should ...
      % go from one to nvar*nvar since the contraint is a big vector.
      Aeq(1,element_restricted) = 1; %zero-impact restriction(s)
end
% NOTICE. For some reasons which I am not fully understannding I can only
% impose independent restrictions of the first column of B, i.e.
% element_restricted can go from 1 to 3. If it is 4, 5, or 6 then the
% second column of B is zero. If it is 7, 8, or 9 then the third column of B is
% zero. My intuition is that it may be related to the long run
% restrictions. The zeros in Xi*B may imply that setting a zero in a specific element of
% some column of B implies that all the column will be zero for
% cointegration. In any case, I need more thoughts on it.

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

%Minimization
B_zero = rand(nvar);
if r > 1
      %B_opt = fmincon(obj, B_zero,[],[],Aeq,Beq,[],[],[],options);
      B_opt = fmincon(obj, B_zero,[],[],Aeq,Beq,[],[],@(B) constraint_long(B,Xi,r,nvar),options);
      %fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
else
      B_opt = fmincon(obj, B_zero,[],[],[],[],[],[],@(B) constraint_long(B,Xi,r,nvar),options);
end
%Temporary tools to visualize the correctness of restrictions
obj_min = obj(B_opt)
obj_zero = obj(B_zero)
B_opt
long_restr = Xi*B_opt

%Setting the A matrix to obtain a SVAR functional form as follows
% y(t) = A1*y(t-1) + ... + Ap*y(t-p) + B*eps(t)
% Again, I am followinf Lutkepohl (2005)
A = zeros(nvar,nvar*(nlags+1));
A(:,1:nvar) = alp*bet' + eye(nvar) + Gam(:,1:nvar);
if nlags >= 2
      for i_lags = 1:nlags-1
            A(:,(i_lags*nvar)+1:(i_lags+1)*nvar) = ...
                  Gam(:,(i_lags*nvar)+1:(i_lags+1)*nvar) ...
                  - Gam(:,((i_lags-1)*nvar)+1:i_lags*nvar);
      end
end
A(:,(nlags*nvar)+1:nvar*(nlags+1)) = - Gam(:,((nlags-1)*nvar)+1:nlags*nvar);

end