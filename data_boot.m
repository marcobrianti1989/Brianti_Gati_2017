function [dataset_boot] = data_boot(B, nburn, res, nsimul)
% Inputs:
% B: Beta coefficient matrix from VAR (nvar*nlags + 1 x nvar)
% nburn = number of burn-in periods
% nlag = number of lags in VAR
% res = residuals from original VAR (Tx1)
% nsimul = number of simulations, i.e. how many bootstrapped datasets to
% generate

T = size(res,1); % time periods
nvar = size(B,2);
nlag = (size(B,1)-1)/nvar;
const = 1;
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
dataset_boot = dataset_boot(nburn+1:end,:,:);

end
