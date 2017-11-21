function [s, obj_opt] = get_structural_shocks(A,gamma,resid)
    
    % A     = the impact matrix recovered through some strategy
    % resid = the seris of reduced form errors from a RF VAR
    
    % Explanation:
    % AA' = sigma (VC matrix)
    % res = (T, nvar) , so:
    % Ass'A' = i'i
    % <=> As = i'
    % => S = A^(-1)i'.
    
    %Set gamma equal to zero if you are doing a plain identification
    %strategy such as short run restrictions, long run. If you are doing a
    %partial identification strategy in the spirit of BarskyandSims then
    %you need to implement the minimization problem and input gamma_opt
    
    T = size(resid,1);
    
    if gamma == 0
        s = A\resid';
    else
        for i_t = 1:T
            resid_t = resid(i_t,:);
            resid_t = resid_t';
            unknown_zero = zeros(size(gamma,2),1);
            objective_zero = obj_get_structural_shocks(A,gamma,resid_t,unknown_zero);
            objective_fun = @(unknown) obj_get_structural_shocks(A,gamma,resid_t,unknown);
            
            %Optimization Parameters
%             options  = optimset('fmincon');
%             options  = optimset(options, 'TolFun', 1e-5, 'display', 'none', 'MaxIter', 1000 );    
            options  = optimoptions('fmincon', 'TolFun', 1e-9, 'display', 'none', 'MaxIter', 1000, ...
                'StepTolerance',1.0000e-15, 'ConstraintTolerance', 1.0000e-09, 'Algorithm', 'sqp' );    
            s(:,i_t) = fmincon(objective_fun, unknown_zero,[],[],[],[],[],[],[],options);
%           s(:,i_t) = fsolve(objective_fun,unknown_zero);
            obj_opt(i_t) = obj_get_structural_shocks(A,gamma,resid_t,s(:,i_t));
        end 
        err = sum(sum(obj_opt.^2));
    end
    
end
    
    
    