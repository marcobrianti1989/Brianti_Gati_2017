function [dataset_boot] = data_boot(B,nburn,res,nsimul,which_correction,q)
% Inputs:
% B: Beta coefficient matrix from VAR (nvar*nlags + 1 x nvar)
% nburn = number of burn-in periods
% nlag = number of lags in VAR
% res = residuals from original VAR (Tx1)
% nsimul = number of simulations, i.e. how many bootstrapped datasets to
% generate
% q = for block drawing, the length of each block

T = size(res,1); % time periods
nvar = size(B,2);
nlag = (size(B,1)-1)/nvar;
const = 1;

switch which_correction
    case 'none'
        dataset_boot = zeros(T+nburn,nvar,nsimul);
        
        for i_repeat = 1:nsimul
            reg_new = zeros(1,nlag*nvar);
            for i_boot = 1:T+nburn
                res_boot = res(randperm(T,1),:);
                yhat = [const, reg_new]*B;
                ystar = yhat + res_boot;
                reg_new = [ystar reg_new(1:end-nvar)];
                dataset_boot(i_boot,:,i_repeat) = ystar;
            end
        end
    case 'blocks'
        block_starts = 1:q:T-mod(T,q); % vector of initial indexes of each block
        k = (T+nburn - mod(T+nburn,q)) ; % how long the sample will be once you discard the leftover
        dataset_boot = zeros(k,nvar,nsimul);

        for i_repeat = 1:nsimul
            reg_new = zeros(1,nlag*nvar);
            for i_boot = 1:q:k % cycles thru T in steps of q
                draw = datasample(block_starts,1); % draw one index from the vector of initial indices...
                for j=1:q % ... and for that index, go thru all shocks in that block
                    res_boot = res(draw+j-1,:); 
                    yhat = [const, reg_new]*B;
                    ystar = yhat + res_boot;
                    reg_new = [ystar reg_new(1:end-nvar)];
                    dataset_boot(i_boot+j-1,:,i_repeat) = ystar;
                end
            end
        end
    otherwise
        error('Correction needs to be either "none", or "blocks".')
        
end

dataset_boot = dataset_boot(nburn+1:end,:,:);

end
