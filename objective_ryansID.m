function [obj_FEV, IRFs_all, FEV_news, FEV_IT] = objective_ryansID(which_variable,H,bet,A,gam)

% which_variable: It select the variable that the shock is maximizing (TFP)
% H: cumulated horizon of the mazimization
% beta is the reduced-form beta and is ((nvar*nlags),nvar) dimension matrix
% A: random impact matrix coming from whatever identification, here we use
% one from cholesky is (nvar*nvar). Importantly, using TFP as a first
% variable with cholesky is automatically imposing constraint (9.5) ->
% A0(1,j) = 0 for all j > 1.
% gamma: is a (nvar x 2) matrix (which is vectorized to be (2*nvar,1) ) which implicitly should select the shocks we are
% interested in. In the main file you need gamma = D(which_shocks) where D =
% eye(nvar). gamma is a selection matrix: it is basically doing the "which_shock"

[nvarlags, nvar]     = size(bet);
nlags                = nvarlags/nvar;
IRFs_all             = zeros(nvar,H);
nshocks              = nvar;

for i_shock = 1:nshocks
    
      shocks              = zeros(nvar,1);
      shocks(i_shock,1)   = 1;      
      % Initialize:
      IRFs_all(:,1,i_shock)   = A*shocks;
      F                       = [IRFs_all(:,1,i_shock)' zeros(1,(nlags-1)*nvar)];
      % Generate IRFs
      for k=2:H
            IRFs_all(:,k,i_shock)    = F*bet;
            F                        = [IRFs_all(:,k,i_shock)' F(1:end-nvar)];
      end
end

DEN = sum(sum(IRFs_all(which_variable,:,:).^2)); %independent from gamma
 
% 1) objective for news shock
% Initialize:
obj_IRFs1(:,1)            = A*gam(1:nvar,1); % gamma3, i.e. corresponding to news
F                         = [obj_IRFs1(:,1)' zeros(1,(nlags-1)*nvar)];
% Generate IRFs
for k=2:H
    obj_IRFs1(:,k)        = F*bet;
    F                     = [obj_IRFs1(:,k)' F(1:end-nvar)];
end

% 2) objective for IT shock
% Initialize:
obj_IRFs2(:,1)            = A*gam(nvar+1:end,1); % gamma4, i.e. corresponding to IT
F                         = [obj_IRFs2(:,1)' zeros(1,(nlags-1)*nvar)];
% Generate IRFs
for k=2:H
    obj_IRFs2(:,k)        = F*bet;
    F                     = [obj_IRFs2(:,k)' F(1:end-nvar)];
end

obj_FEV  = sum(obj_IRFs1(which_variable,:).^2) + sum(obj_IRFs2(which_variable,:).^2);
obj_FEV  = obj_FEV/DEN; %normalization
FEV_news = sum(obj_IRFs1(which_variable,:).^2)/DEN;
FEV_IT   = sum(obj_IRFs2(which_variable,:).^2)/DEN;


obj_FEV = - obj_FEV; %BE CAREFUL! IT IS ALREADY NEGATIVE FOR FMINCON!!!