function [alph_hat, bet_hat,Gam_hat,sigma,obj] = RRR(dataset, nlags, r)

% INPUTS:
% dataset is (T, nvar)
% nlags sets the number of lags (scalar)
% r implictly decides how many cointegrating vectors we have. (scalar)
% Number of cointegrating vectors is (r)
% In this case we have (nvar-r) permanent shocks and r temporary shocks

% OUTPUT:
% Rearrange data to run a reduced rank regression, outputting:
% bet_hat is the cointegrating vector and is (nvar, r)
% alph_hat is the loading vector and is (nvar, r)
% Gam_hat is the betas of the lagged differentiated terms. Similar to Reduced
% form VAR and is (nvar, nvar*nlags) - No constant for now!
% sigma the variance covariance matrix of res (nvar, nvar)

% Comment. With T' I mean T properly adjusted for nlags.

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
% an anonymous document saved as ReducedRankRegression.pdf

% y = alp*bet'*x + Gam*z + U
% where y is (nvar,T'), alp is (nvar,r), bet is (nvar,r), x is (nvar,T'),
% Gam is (nvar,nvar*nlags), z is (nvar*nlags,T') and U is (nvar,T')
x      = dataset(nlags+1:end-1,:)'; 
% I am removing the last obs from Y since it is x which implies that
% Y(T) will never be used. 
z      = reg'; %dT(t-j) first diff of the right-hand side
y      = diff_dataset(nlags+1:end,:)'; %left-hand side - Should be equal to data_1
% dY should be equal to data_1. Otherwise there is something wrong

%Following the algorithm
Syx = (y*x')/T;
Sxy = (x*y')/T;
Sxx = (x*x')/T;
Syy = (y*y')/T;
Syz = (y*z')/T;
Szy = (z*y')/T;
Szz = (z*z')/T;
Sxz = (x*z')/T;
Szx = (z*x')/T;

Syxz = Syx - Syz*((Szz)^(-1))*Szx;
Sxyz = Sxy - Sxz*((Szz)^(-1))*Szy;
Syyz = Syy - Syz*((Szz)^(-1))*Szy;
Sxxz = Sxx - Sxz*((Szz)^(-1))*Szx;

A = Sxyz*((Syyz)^(-1))*(Syxz);
B = Sxxz;

eigv = eig(A,B);
sorted_eig = sort(eigv,'descend');

%Check if the generalize eigenvalue problem is correct
[V, D] = eig(A,B); %V is the eigenvectors matrix and D is the eignvalue diagonal matrix
%[V,D] = eig(A,B) produces a diagonal matrix D of generalized
%eigenvalues and a full matrix V whose columns are the corresponding
%eigenvectors so that A*V = B*V*D.
check = A*V - B*V*D;
check = sum(sum(abs(check)));
if check > 10^(-12)
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
alph_hat = Syxz*bet_hat;

%Getting Gamma_hat
Gam_hat = (y - alph_hat*bet_hat'*x)*z'*(z*z')^(-1); %(nvar,nvar*nlags)

%Getting residuals and Variance-Covariance Matrix 
sigma = Syyz - Syxz*bet_hat*(((bet_hat')*Sxxz*bet_hat)^(-1))*(bet_hat')*Sxyz;

%y = alp*bet'*x + Gam*z + res
res = alph_hat*bet_hat'*x + Gam_hat*z - y;
obj = sum(sum(res.^2));

end


