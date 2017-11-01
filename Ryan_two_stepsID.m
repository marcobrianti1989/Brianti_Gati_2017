function [impact, FEV_opt, IRFs, gam3_opt, FEV_news, FEV_IT] ...
        = Ryan_two_stepsID(which_variable,which_shocks,H,B,A,q)
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
% gamma0 = gamma0(:);
gam3_zero = gamma0(:,1);
gam4_zero = gamma0(:,2);

%First Step - Identifying gam3 - Similar to B&S(2012) with long-restriction
%over relative price

% Setting the objective function of Step1
obj = @(gam3) objective_barskysims(which_variable,H,B,A,gam3);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

% LR constraint that news has no effect on P(IT)/P in the LR.
G = eye(nvar);
for i = 1:nlags
    G = G - B((i-1)*nvar+1:i*nvar,:);
end
P = inv(G); % LR effect overall
R = P*A; %Long-run IR matrix a la Blanchard and Quah

%Constraint that News and IT have no contemporaneous effect on TFP
Me3       = 1; %  Me = no. of equality constraints
Beq3      = zeros(Me3,1); % Beq is (Me x 1) where
Aeq3      = zeros(Me3,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq3(1,1) = 1; %zero-impact of news on TFP
% Aeq3(2,q) = R(q,3); % LR restriction - wrong!

% [gam3_opt] = fmincon(obj, gam3_zero,[],[],Aeq3,Beq3,[],[],@(gam3) constraint_barskysims11(gam3),options);
[gam3_opt] = fmincon(obj, gam3_zero,[],[],Aeq3,Beq3,[],[],@(gam3) constraint_ryan(gam3,R,q),options);

%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
[FEV3_opt, ~, ~]     = objective_barskysims(which_variable,H,B,A,gam3_opt);
FEV3_opt     = - FEV3_opt;

if FEV3_opt > 1 || (gam3_opt'*gam3_opt - 1)^2 > 10^(-14) || gam3_opt(1)^2 > 10^(-14)
    error('The problem is not consistent with the constraints.')
end

%Second Step - Identifying gam4 - Max remaining FEV_TFP using IT shock

% Setting the objective function of Step2
obj = @(gam4) objective_RyanID_twosteps(which_variable,H,B,A,gam3_opt,gam4);

%Constraint that News and IT have no contemporaneous effect on TFP
Me4       = 1; %  Me = no. of equality constraints
Beq4      = zeros(Me4,1); % Beq is (Me x 1) where
Aeq4      = zeros(Me4,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq4(1,1) = 1; %zero-impact restrictions of news on TFP

[gam4_opt] = fmincon(obj, gam4_zero,[],[],Aeq4,Beq4,[],[],@(gam4) constraint_orthogonality_ryan([gam3_opt gam4]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
[FEV_opt,IRFs,FEV_news,FEV_IT] = objective_RyanID_twosteps(which_variable,H,B,A,gam3_opt,gam4_opt);
FEV_opt     = - FEV_opt;

if FEV_opt > 1 || (gam4_opt'*gam4_opt - 1)^2 > 10^(-14) || gam4_opt(1)^2 > 10^(-14)
    error('The problem is not consistent with the constraints.')
end

gam = [gam3_opt gam4_opt];
impact = A*gam;

end




