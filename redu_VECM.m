function [alph_hat,bet_hat,Pi,Gam_hat,U,sigma] = redu_VECM(dataset, nlags, r)

% INPUTS:
% dataset is (T, nvar)
% nlags sets the number of lags (scalar)
% r implictly decides how many cointegrating vectors we have. (scalar)
% Number of cointegrating vectors is (r)
% In this case we have (nvar-r) permanent shocks and r temporary shocks

% OUTPUT:
% Rearrange data to run a reduced form VECM, outputting:
% bet_hat is the cointegrating vector and is (nvar, r)
% alph_hat is the loading vector and is (nvar, r)
% Pi, the correction matrix. (nvar, nvar)
% Gam_hat is the betas of the lagged differentiated terms. Similar to Reduced
% form VAR and is (nvar, nvar*nlags) - No constant for now!
% res = residuals from VECM estimation. (nvar,T')
% sigma the variance covariance matrix of res (nvar, nvar)

% Comment. With T' I mean T properly adjusted for nlags. (T'= T-nlags-1)

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

% Defining objects - For now without constant. Notations is consistent with
% Lutkepohl (2005) - EUI Working Paper

% See Equation (3.1). The equation is the following
% dY = alp*bet'*Y + Gam*X + U
% where dY is (nvar,T'), alp is (nvar,r), bet is (nvar,r), Y is (nvar,T'),
% Gam is (nvar,nvar*nlags), X is (nvar*nlags,T') and U is (nvar,T')
Ylag       = dataset(nlags+1:end-1,:)'; %Y(t-1) level of the right-hand side.
% I am removing the last obs from Y since it is Y(t-1) which implies that
% Y(T) will never be used. 
X       = reg'; %dT(t-j) first diff of the right-hand side
dY      = diff_dataset(nlags+1:end,:)'; %left-hand side - Should be equal to data_1
% dY should be equal to data_1. Otherwise there is something wrong
M       = eye(T-nlags-1,T-nlags-1) - X'*(X*X')^(-1)*X; %(T',T')
S00     = (dY*M*dY')/T; %(nvar,navr)
S01     = (dY*M*Ylag')/T;  %(nvar,navr)
S11     = (Ylag*M*Ylag')/T;   %(nvar,navr)

% e = eig(A,B) returns a column vector containing the ...
% generalized eigenvalues of square matrices A and B.
% det(lambda*B - A)
A = S01'*(S00^(-1))*S01;
B = S11;
eigv = eig(A,B);
sorted_eig = sort(eigv,'descend');

%Check if the generalize eigenvalue problem is correct
[V, D] = eig(A,B); %V is the eigenvectors matrix and D is the eignvalue diagonal matrix
%[V,D] = eig(A,B) produces a diagonal matrix D of generalized
%eigenvalues and a full matrix V whose columns are the corresponding
%eigenvectors so that A*V = B*V*D.
check = A*V - B*V*D;
check = sum(sum(abs(check)));
if check > 10^(-14)
      warning('Eigenvalues and eigenvectors are not solving the generalized eigenvalue problem.')
end

%Resorting the eigenvectors using the ordering of descedent eignvalues
bet_full = zeros(size(V,1),size(V,2));
for i_eig = 1:nvar
      idx = find(abs(sorted_eig(i_eig) - diag(D))<=eps);
      bet_full(:,i_eig) = V(:,idx);
end

%Getting alpha_hat and beta_hat
bet_hat = bet_full(:,1:r);
% alph_hat = dY*M*Ylag'*bet_hat*(bet_hat'*Ylag*M*Ylag'*bet_hat)^(-1);
alph_hat = S01*bet_hat; % Prior: this one is correct, line 85 is wrong. Check!
Pi = alph_hat*bet_hat';

%Getting Gamma_hat
Gam_hat = (dY - alph_hat*bet_hat'*Ylag)*X'*(X*X')^(-1); %(nvar,nvar*nlags)

%Getting residuals and Variance-Covariance Matrix
U = dY*M - alph_hat*bet_hat'*Ylag*M; 
sigma = U*U'/(T-nlags-1);

U2 = - alph_hat*bet_hat'*Ylag - Gam_hat*X + dY;
obj = sum(sum(U2.^2));

end