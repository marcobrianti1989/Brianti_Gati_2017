function A = barskysims(which_variable,which_shock,H,B,A)
%which_variable: variable we want to max the response for a specific shock
%which_shock: the shock is maximizing the FEV of which_variable
%H: horizon of the maximization
%B: B ols with constant from the reduform of the VAR is (nvar*nlag+1,nvar)
%A: initial identified matrix. For B&S we use standar chol with TFP first

B                    = B(2:end,:);
[nvarlags, nvar]     = size(B);
nlags                = nvarlags/nvar;
    
%Barsky and Sims identification strategy
D = eye(nvar);
gamma0 = D(:,which_shock);

% dbstop in objective_barskysims at 21
[FEV] = objective_barskysims(which_variable,H,B,A,gamma0);
        
obj = @(gam) objective_barskysims(which_variable,H,B,A,gam);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'iter');

%Constraint 10 from Barsky and Sims, News have no contemporaneous effect on TFP
Beq      = zeros(nvar,1);
Aeq      = zeros(nvar,nvar);
Aeq(1,1) = 1;

[gamma_opt] = fmincon(obj, gamma0,[],[],Aeq,Beq,[],[],@(gam) constraint_barskysims11(gam),options);
             %fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
FEV_opt     = objective_barskysims(which_variable,H,B,A,gamma_opt);
FEV_opt     = - FEV_opt;

if FEV_opt > 1 || (gamma_opt'*gamma_opt - 1)^2 > 10^(-14) || gamma_opt(1)^2 > 10^(-14) 
    error('The problem is not consistent with the constraints.')
end

A(:,which_shock) = A*gamma_opt;











