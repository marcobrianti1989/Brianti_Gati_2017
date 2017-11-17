function [B,res,sigma] = reduform_var(dataset, nlags)

% Rearrange data to run a reduced-form VAR, outputting B, the reduced fprm
% betas.
% res = residuals from VAR estimation
% sigma = VC matrix of RF VAR

T = size(dataset,1); % time periods
nvar = size(dataset,2);

%Taking the lag on the dataset
for ilag = 1:nlags+1
    eval(['data_', num2str(ilag), ' = dataset(2+nlags-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(T-nlags,nlags*nvar);
for ilag = 2:nlags+1
    for ivar = 1:nvar
        eval(['reg(:,ivar+size(dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
    end
end

% Add constant
reg = [ones(T-nlags,1) reg]; % (T-nlags, 1+nvar*nlags)

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)\(reg'*data_1);

%Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

%Compute the variance-covariance matrix
sigma = res'*res/(T-nlags);