function [Pi,B,res,sigma,bet,alph] = VECM(dataset, nlags, constant, r)


% INPUTS: 
% dataset is (T, nvar)
% nlags sets the number of lags (scalar)
% constant equal to 1 the VECM hat the constant, different from one no
% constant. (scalar)
% r decides how many cointegrating vectors we have. (scalar)

% OUTPUT:
% Rearrange data to run a VECM, outputting Pi, the correction matrix. (nvar, nvar)
% B is the betas of the lagged differentiated terms. Similar to Reduced
% form VAR. B is ((nvar*nlags) + 1, nvar)
% res = residuals from VECM estimation. (T - nlags - 1, nvar)
% sigma the variance covariance matrix of res (nvar, nvar)
% bet is the cointegrating vector and is (nvar, r)
% alph is the loading vector and is (nvar, r)

T               = size(dataset,1); % time periods
nvar            = size(dataset,2); % number of variables
diff_dataset    = diff(dataset); % taking the first difference of the whole dataset
 
% Taking the lag on the dataset
for ilag = 1:nlags+1
      eval(['data_', num2str(ilag), ' = diff_dataset(2+nlags-ilag:end+1-ilag,:);'])
end

% Building the matrix of regressor (which are the lags of the left hand side)
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

% OLS in the Reduced Form Vector Autoregression
B = (reg'*reg)\(reg'*data_1);

% Evaluate the residuals (res = yhat - y)
res = data_1 - reg*B;

% Compute the variance-covariance matrix
sigma = res'*res/(T-nlags);

% Compute the correction matrix
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


r = 5;
% Set random vector beta in order to split Pi in alpha times beta'.
%bet = [ones(1,r); ones(nvar-1,r)];
%Set the objective function to min over alpha
obj_alph = @(alph_bet_x) objective_alph(Pi,r,alph_bet_x);
%Initial value
alph_bet_zero = ones(nvar,2*r);

%Optimization Parameters
options  = optimset('fmincon');
options  = optimset(options, 'TolFun', 1e-9, 'display', 'none');
%Optimization
[alph_bet_opt] = fmincon(obj_alph, alph_bet_zero,[],[],[],[],[],[],[],options);
%fmincon(FUN,X0,A,B,Aeq,Beq,LB,UB,NONLCON,OPTIONS)
%alph_bet_opt = fzero(obj_alph,alph_bet_zero);


alph = alph_bet_opt(:,r);
bet  = alph_bet_opt(:,r+1:end);

obj_alph(alph_bet_opt)

Pi - alph*bet'





end