function [impact, FEV_opt, IRFs, gamma_opt, FEV_news, FEV_IT] = ryansID(which_variable,which_shocks,H,B,A,q)
%which_variable: variable we want to max the response for a specific shock (TFP)
%which_shocks: the shock is maximizing the FEV of which_variable (--> news and IT shocks)
%H: horizon of the maximization
%B: B ols with constant from the reduform of the VAR is (nvar*nlag+1,nvar)
%A: initial identified matrix. Also here we use standar chol with TFP first
% q: position of rel. price of IT

B                    = B(2:end,:);
[nvarlags, nvar]     = size(B);
nlags                = nvarlags/nvar;
nshocks              = length(which_shocks);

%Barsky and Sims identification strategy
D = eye(nvar);
gamma0 = D(:,which_shocks);
gamma0 = gamma0(:);

% dbstop in objective_ryansID at 21
[FEV, IRFs] = objective_ryansID(which_variable,H,B,A,gamma0);

obj = @(gam) objective_ryansID(which_variable,H,B,A,gam);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

% LR constraint that news has no effect on P(IT)/P in the LR.
G = eye(nvar);
for i = 1:nlags
    G = G - B((i-1)*nvar+1:i*nvar,:);
end
P = inv(G); % LR effect overall

%Constraint that News and IT have no contemporaneous effect on TFP
Me = 3; %  Me = no. of equality constraints
Beq      = zeros(Me,1); % Beq is (Me x 1) where
Aeq      = zeros(Me,nshocks*nvar); % Aeq is (Me x (nshocks*nvar))
Aeq(1,1) = 1;
Aeq(2,7) = 1;
Aeq(3,q) = P(q,3)*A(q,3); % LR restriction
% Aeq(3,end) = 1; % zero impact restriction of IT on rel. prices

[gamma_opt] = fmincon(obj, gamma0,[],[],Aeq,Beq,[],[],@(gam) constraint_orthogonality_ryan(gam),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
[FEV_opt, ~, FEV_news, FEV_IT]     = objective_ryansID(which_variable,H,B,A,gamma_opt);
FEV_opt     = - FEV_opt;

if FEV_opt > 1 || (gamma_opt'*gamma_opt - 1)^2 > 10^(-14) || gamma_opt(1)^2 > 10^(-14)
    error('The problem is not consistent with the constraints.')
end

gamma_opt = [gamma_opt(1:nvar,1), gamma_opt(nvar+1:end,1)];
impact = A*gamma_opt;
