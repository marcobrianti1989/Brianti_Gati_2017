function [Pi,B,res,sigma] = VECM(dataset, nlags, constant)

% Rearrange data to run a VECM, outputting Pi, the correction matrix.
% reg the betas of the lagged terms. Similar to Reduced form VAR
% res = residuals from VAR estimation
% sigma the variance covariance matrix.
% constant equal to 1 the VECM hat the constant, different from one no
% constant.

T = size(dataset,1); % time periods
nvar = size(dataset,2);
diff_dataset = diff(dataset);

%Taking the lag on the dataset
for ilag = 1:nlags+1
      eval(['data_', num2str(ilag), ' = diff_dataset(2+nlags-ilag:end+1-ilag,:);'])
end

%Building the matrix of regressor (which are the lags of the left hand side)
reg = zeros(T-nlags-1,nlags*nvar);
for ilag = 2:nlags+1
      for ivar = 1:nvar
            eval(['reg(:,ivar+size(diff_dataset,2)*(ilag-2)) = data_',num2str(ilag),'(:,ivar);'])
      end
end

% Add/donotadd constant and correction term
if constant == 1
      reg = [ones(T-nlags-1,1) dataset(3:end,:) reg];
else
      reg = [dataset(3:end,:) reg];
end

%OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)\(reg'*data_1);

%Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

%Compute the variance-covariance matrix
sigma = res'*res/(T-nlags);

%Compute the correction matrix
if constant == 1
      Pi   = B(2:nvar+1,:);
      B    = [B(1,:);   B(nvar+2:end,:)];
else
      Pi   = B(1:nvar,:);
      B    = [B(nvar+1:end,:)];
end

if det(Pi) >= 10^(10-8)
      warning('Determinant of cointegration matrix is not numerically zero.')
end

%Set vector beta in order to split Pi in alpha times beta'.
y_var = dataset(:,1);
reg_x = dataset(:,2:end);

beta = (reg_x'*reg_x)\(reg_x'*y_var);
beta = -beta;
beta = [1; beta];



end