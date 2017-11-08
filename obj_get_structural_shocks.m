function objective_fun = obj_get_structural_shocks(A,gamma,resid,unknowns)
    %For the case with partial identification similar to Barsky&uncleSims
    %we do not have D such that ADD'A = sigma then we basically solve with
    %fmin con the following problem: A*gamma*structural_t = resid_t
    %resid is (1,nvar), A is (nvar,nvar), gamma is (nvar,size(gamma,2)) in our case,
    %and structural (which is unknown variable) is (size(gamma,2), 1)
    resid_hat = A*gamma*unknowns;
    objective_fun = sum((resid_hat - resid).^2);
    
end