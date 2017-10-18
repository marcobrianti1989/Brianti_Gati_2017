function [] = objective_barskysims(which_shock,which_variable,H,beta,A,D) 
% Equation 7 from Barsky and Sims (2011)
% which_shock: It select the shock whose effect wants to be maximized
% which_variable: It select the variable that the shock is maximizing
% H: cumulated horizon of the mazimization
% beta is the reduced-form beta and is ((nvar*nlags)+1,nvar) dimension matrix
% A: random impact matrix coming from whatever identification, here we use
% one from cholesky is (nvar*nvar)

beta              = beta(2:end,:);
[nvarlags, nvar]  = size(beta);
nlags             = nvarlags/nvar;

e_s                   = zeros(nvar,1);
e_v                   = zeros(nvar,1);
e_s(which_shock)      = 1;
e_s(which_variable)   = 1;

gamma = D*e_s






end