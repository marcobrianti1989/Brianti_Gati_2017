function [B_opt, Xi, A, diff_BBp_sigma] = structural_VECM(alp,bet,Gam,res,sigma,nlags,r)

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

% NOTICE:
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

%We need to add 1/2*nvar(nvar - 1) - r*(nvar - r) zero impact restrictions
%In particular,
%1. r*(r - 1)/2 must be add on temporary shocks
%2. (nvar - r)*(nvar - r - 1)/2 must be add on permanent shocks
%Me1: number of impact restrictions for temporary shocks
%Me2: number of impact restrictions for permanent shocks
if r > 1 && nvar > r + 1 %We need restrictions for both permanent and temporary
      Me1                                    = r*(r-1)/2; %  Me = no. of equality constraints.
      Me2                                    = (nvar-r)*(nvar-r-1)/2;
      Me                                     = Me1 + Me2; %  Me = no. of equality constraints.
      Beq                                    = zeros(Me,1); % Beq is (Me x 1)
      Aeq                                    = zeros(Me,nvar*nvar); % Aeq is (Me x (nvar*nvar)) - it is B vectorized
      Aeq(1:Me2,1:(nvar-r)*(nvar-r-1)/2)     = eye(Me2); %zero-impact restriction(s) permanent shocks
      Aeq(Me2+1:Me,end-Me1+1:end)            = eye(Me1); %zero-impact restriction(s) temporary shocks
elseif r > 1 && nvar <= r + 1 %We need restrictions only for temporary
      Me1                                    = r*(r-1)/2; %  Me = no. of equality constraints.
      Me2                                    = 0;
      Me                                     = Me1 + Me2;
      Beq                                    = zeros(Me,1); % Beq is (Me x 1)
      Aeq                                    = zeros(Me,nvar*nvar); % Aeq is (Me x (nvar*nvar)) - it is B vectorized
      Aeq(Me2+1:Me,end-Me1+1:end)            = eye(Me); %zero-impact restriction(s)
elseif r == 1 && nvar > r + 1 %We need restrictions only for permanent
      Me1                                    = 0;
      Me2                                    = (nvar-r)*(nvar-r-1)/2;
      Me                                     = Me1 + Me2; %  Me = no. of equality constraints.
      Beq                                    = zeros(Me,1); % Beq is (Me x 1)
      Aeq                                    = zeros(Me,nvar*nvar); % Aeq is (Me x (nvar*nvar)) - it is B vectorized
      regular = 0;
      if regular == 1
            Aeq(1:Me2,1:(nvar-r)*(nvar-r-1)/2)     = eye(Me2); %zero-impact restriction(s) permanent shocks
      else
            disp('In the VECM you are doing something identification-specific')
            Aeq(1,1) = 1;
            Aeq(2,nvar+1) = 1;
            Aeq(3,nvar*2+1) = 1;
      end
      
      %Optimization Parameters
      options  = optimset('fmincon');
      options  = optimset(options,'TolFun', 1e-19);
      %'TolFun', 1e-19, 'FinDiffRelStep', 1
      warning off
      %Minimization
      B_zero = randn(nvar,nvar);
      if r == 1 && nvar <= r + 1 %We do not need impact restrictions
            B_opt = fmincon(obj, B_zero,[],[],[],[],[],[],@(B) constraint_long(B,Xi,r,nvar),options);
      else %We need some sort of impact restrictions
            B_opt  = fmincon(obj, B_zero,[],[],Aeq,Beq,[],[],@(B) ...
                  constraint_long(B,Xi,r,nvar),options);
            %fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
      end
      
      %Temporary tools to visualize the correctness of restrictions
      diff_BBp_sigma = sum(sum((B_opt*B_opt' - sigma).^2));
      long_restr = Xi*B_opt;
      B_opt;
      warning on
      
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