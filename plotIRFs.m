function plotIRFs(A,A_boot,B,B_boot,h,which_shock, names, varnames, sig)
% h = IR horizon
% which_shock should be input as a vector e.g. [2 3] meaning the second and
% third shock.
% sig = CI significance level (enter as 0.9 for example)
% names = cell vector of shock names

nvar = size(A,1);
nlag = (size(B,1)-1)/nvar;
nshocks = size(which_shock,2);
periods = 1:h;
nsimul = size(A_boot,3);
perc_up = ceil(nsimul*sig); % the upper percentile of bootstrapped responses for CI
perc_low = floor(nsimul*(1-sig)); % the lower percentile of bootstrapped responses for CI

% Remove constant
B = B(2:end,:);
B_boot = B_boot(2:end,:,:);

for i_shock=1:nshocks
    shocks = zeros(nvar,1);
    shocks(which_shock(i_shock),1) = 1;
    name = names{i_shock};
    
    % Preallocate IRs
    IRFs = zeros(nvar,h,nshocks);
    IRFs_boot = zeros(nvar,h,nshocks,nsimul);
    IRFs_boot_sorted = zeros(nvar,h,nshocks,nsimul);
    ub = zeros(nvar,h,nshocks);
    lb = zeros(nvar,h,nshocks);
    
    % Initialize:
    IRFs(:,1,i_shock) = A*shocks;
    F = [IRFs(:,1,i_shock)' zeros(1,(nlag-1)*nvar)];
    % Generate IRFs
    for k=2:h
        IRFs(:,k,i_shock) = F*B;
        F = [IRFs(:,k,i_shock)' F(1:end-nvar)];
    end
    
    % Redo for bootstrapped samples: Here we need an extra loop over nsimul
    for i_sim=1:nsimul
%       Initialize:
        IRFs_boot(:,1,i_shock,i_sim) = A*shocks;
        F = [IRFs_boot(:,1,i_shock,i_sim)' zeros(1,(nlag-1)*nvar)];
        % Generate IRFs
        for k=2:h
            IRFs_boot(:,k,i_shock,i_sim) = F*B;
            F = [IRFs_boot(:,k,i_shock,i_sim)' F(1:end-nvar)];
        end
    end
    
    % Sort bootstrap IRFs and set lower and upper bounds
    for i_shocks = 1:nshocks
        IRFs_boot_sorted(:,:,i_shocks,:) = sort(IRFs_boot_sorted(:,:,i_shocks,:),4);
        ub(:,:,i_shocks) = IRFs_boot_sorted(:,:,i_shocks,perc_up);
        lb(:,:,i_shocks) = IRFs_boot_sorted(:,:,i_shocks,perc_low);
    end
    % TO DO: add ub and lb on graphs!
    
    % Draw pretty pictures
    figure(i_shock)
    for i_var=1:nvar
        varname = varnames{i_var};
        subplot(nvar,1,i_var)
        plot(periods,IRFs(i_var,:,i_shock),'linewidth',2)
        title([name, ' on ' , varname])
        grid on
        
    end
    
end