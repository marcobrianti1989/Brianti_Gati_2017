function [impact, impact_IT_opt, gam_opt]= just_IT_ID(which_variable,which_shock,A)

%which_variable: variable we want to max the response for a specific shock (TFP)
%which_shock: the shock is maximizing the FEV of which_variable (--> IT shock)
%A: initial identified matrix. Also here we use standar chol with TFP first

nvar = size(A,2); 

%Barsky and Sims identification strategy
D = eye(nvar);
gam_zero = D(:,which_shock); %IT shock impact vector (initial value)

%% Identifying gam4 - Max impact FEV_IT s.t. impact on TFP = 0

% Setting the objective function of Step2
obj = @(gam) objective_just_IT(which_variable,A,gam);

%Constraint that News and IT have no contemporaneous effect on TFP
Me       = 1; %  Me = no. of equality constraints
Beq      = zeros(Me,1); % Beq is (Me x 1) where
Aeq      = zeros(Me,1*nvar); % Aeq is (Me x (nshock*nvar)) - nshock is 1 at this step
Aeq(1,1) = 1; %zero-impact restrictions of IT on TFP

gam_opt = fmincon(obj, gam_zero,[],[],Aeq,Beq,[],[],@(gam) constraint_orthogonality_ryan(gam));
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
impact_IT_opt = objective_just_IT(which_variable,A,gam_opt);
impact_IT_opt     = - impact_IT_opt;

if impact_IT_opt > 1 || (gam_opt'*gam_opt - 1)^2 > 10^(-10) || gam_opt(1)^2 > 10^(-12)
    warning('The problem is not consistent with the constraints.')
end

impact = A*gam_opt;

end

