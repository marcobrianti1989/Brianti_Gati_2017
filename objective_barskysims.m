function obj_FEV = objective_barskysims(which_variable,H,beta,A,gam) 
% Equation 7 from Barsky and Sims (2011)
% which_shock: It select the shock whose effect wants to be maximized
% which_variable: It select the variable that the shock is maximizing
% H: cumulated horizon of the mazimization
% beta is the reduced-form beta and is ((nvar*nlags)+1,nvar) dimension matrix
% A: random impact matrix coming from whatever identification, here we use
% one from cholesky is (nvar*nvar). Importantly, using TFP as a first
% variable with cholesky is automatically imposing constraint (9.5) ->
% A0(1,j) = 0 for all j > 1.
% gamma: is a vector which implicitly should select the shock we are
% interested in. In the main file you need gamma = D(which_shock) where D =
% eye(nvar). gamma is a selection vector: it is basically doing the "which_shock"

beta                 = beta(2:end,:);
[nvarlags, nvar]     = size(beta);
nlags                = nvarlags/nvar;
IRFs                 = zeros(nvar,H);
nshocks              = nvar;

for i_shock = 1:nshocks
    
      shocks              = zeros(nvar,1);
      shocks(i_shock,1)   = 1;      
      % Initialize:
      IRFs(:,1,i_shock)   = A*shocks;
      F                   = [IRFs(:,1,i_shock)' zeros(1,(nlags-1)*nvar)];
      % Generate IRFs
      for k=2:H
            IRFs(:,k,i_shock)    = F*beta;
            F                    = [IRFs(:,k,i_shock)' F(1:end-nvar)];
      end
end

DEN = sum(sum(IRFs(which_variable,:,:).^2)); %independent from gamma
  
% Initialize:
obj_IRFs(:,1)   = A*gam;
F                   = [obj_IRFs(:,1)' zeros(1,(nlags-1)*nvar)];
% Generate IRFs
for k=2:H
    obj_IRFs(:,k)    = F*beta;
    F                    = [obj_IRFs(:,k)' F(1:end-nvar)];
end

obj_FEV = sum(obj_IRFs(which_variable,:).^2);
obj_FEV = obj_FEV/DEN; %normalization

obj_FEV = - obj_FEV; %BE CAREFUL! IT IS ALREADY NEGATIVE FOR FMINCON!!!










end