function series = counterfactual(s,imp,B,y0,which_shock_zero)
    % s is the whole matrix of structual shocks (T,nvar)
    % impact is the impact matrix from the SVAR (nvar,nvar)
    % B is the estimated beta OLS from the reduform VAR (1+nvar*nlags,nvar)
    % y0 is the initial value of the interest variable. It is a scalar
    % which_shock zero is a scalar between 1 and nvar and points out which
    % shock should be set to zero
    
    % Technical Values
    [T, nvar] = size(s);
    for it = 1:T
        if s(it,which_shock_zero) < 0
    s(it,which_shock_zero) = 0;
        end
    end
    nlag = (size(B,1)-1)/nvar;
    
    % Initialize:
    series = zeros(T,nvar);
    series(1,:) = y0;
    reg = [ones(1,1) series(1,:) zeros(1,(nlag-1))*nvar];
    % Generate IRFs
    for t= 1:T
        series(t+1,:) = reg*B + (imp*s(t,:)')';
        reg = [ones(1,1) series(t+1,:) reg(1:end-nvar-1)];
    end
    
    
    
    
end
