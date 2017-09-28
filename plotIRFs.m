function plotIRFs(A,A_boot,B,B_boot,h,which_shock, names, varnames)
% h = IR horizon
% which_shock should be input as a vector e.g. [2 3] meaning the second and
% third shock.
% names = cell vector of shock names

nvar = size(A,1);
nlag = (size(B,1)-1)/nvar;
nshocks = size(which_shock,2);
periods = 1:h;

% Remove constant
B = B(2:end,:);
B_boot = B_boot(2:end,:,:);


for i_shock=1:nshocks
    shocks = zeros(nvar,1);
    shocks(which_shock(i_shock),1) = 1;
    name = names{i_shock};
    
    
    % Create IRs
    IRFs = zeros(nvar,h,nshocks);
    
    % Initialize:
    IRFs(:,1,i_shock) = A*shocks;
    F = [IRFs(:,1,i_shock)' zeros(1,(nlag-1)*nvar)];
    for k=2:h
        IRFs(:,k,i_shock) = F*B;
        F = [IRFs(:,k,i_shock)' F(1:end-nvar)];
    end
    
    figure(i_shock)
    for i_var=1:nvar
        varname = varnames{i_var};
        subplot(nvar,1,i_var)
        plot(periods,IRFs(i_var,:,i_shock),'linewidth',2)
        title([name, ' on ' , varname])
        grid on
        
    end
    
end