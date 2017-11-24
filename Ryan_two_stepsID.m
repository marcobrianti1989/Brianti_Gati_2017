function [impact, FEV_opt, IRFs, gam, FEV_news, FEV_IT] ...
        = Ryan_two_stepsID(which_variable,which_shocks,H,LR_hor,B,A,q)
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
gam3_zero = gamma0(:,1); %news shock impact vector (initial value)
gam4_zero = gamma0(:,2); %IT shock impact vector (initial value)

%% First Step - Identifying gam3 - Similar to B&S(2012) with long-restriction
%over relative price

% Setting the objective function of Step1
obj = @(gam3) objective_barskysims(which_variable,H,B,A,gam3); % usual is BS
% obj = @(gam3) objective_IRmax(which_variable,H,B,A,gam3); % try Forni-style IR-max

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');

% Instead of Blanchard Quah, impose that IR(rel_prices, news) = 0 at some
% finite long horizon.
% LR_hor = 8; % 100
sig = 0.9; % not used, but we need to input something not to get error
[IR, ~, ~] = genIRFs(A,0,vertcat(ones(1,size(B,2)), B),0,LR_hor, sig);
% R = zeros(nvar,nvar); 
R = squeeze(IR(:,end,:)); % comment this out if you wanna not impose the LR-restriction

%Constraint that News and IT have no contemporaneous effect on TFP
Me3       = 1; %  Me = no. of equality constraints
Beq3      = zeros(Me3,1); % Beq is (Me x 1) where
Aeq3      = zeros(Me3,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq3(1,1) = 1; %zero-impact of news on TFP

% dbstop constraint_ryan at 12
% [gam3_opt] = fmincon(obj, gam3_zero,[],[],Aeq3,Beq3,[],[],@(gam3) constraint_barskysims11(gam3),options);
[gam3_opt] = fmincon(obj, gam3_zero,[],[],Aeq3,Beq3,[],[],@(gam3) constraint_ryan(gam3,R,q),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)

% Recover FEV from the step 1 minimization
[FEV3_opt, ~, ~]     = objective_barskysims(which_variable,H,B,A,gam3_opt);
FEV3_opt     = - FEV3_opt;

if FEV3_opt > 1 || (gam3_opt'*gam3_opt - 1)^2 > 10^(-10) || gam3_opt(1)^2 > 10^(-12)
    warning('The problem is not consistent with the constraints.')
end

%% Second Step - Identifying gam4 - Max remaining FEV_TFP using IT shock

% Setting the objective function of Step2
obj = @(gam4) objective_RyanID_twosteps(which_variable,H,B,A,gam3_opt,gam4);

%Constraint that News and IT have no contemporaneous effect on TFP
Me4       = 1; %  Me = no. of equality constraints
Beq4      = zeros(Me4,1); % Beq is (Me x 1) where
Aeq4      = zeros(Me4,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq4(1,1) = 1; %zero-impact restrictions of IT on TFP

[gam4_opt] = fmincon(obj, gam4_zero,[],[],Aeq4,Beq4,[],[],@(gam4) constraint_orthogonality_ryan([gam3_opt gam4]),options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
[FEV_opt,IRFs,FEV_news,FEV_IT] = objective_RyanID_twosteps(which_variable,H,B,A,gam3_opt,gam4_opt);
FEV_opt     = - FEV_opt;

if FEV_opt > 1 || (gam4_opt'*gam4_opt - 1)^2 > 10^(-10) || gam4_opt(1)^2 > 10^(-12)
    warning('The problem is not consistent with the constraints.')
end

gam = [gam3_opt gam4_opt];
impact = A*gam;

end




