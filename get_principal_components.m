function pc = get_principal_components(data)
% Get principal components (see Stock & Watson 2002, eq 6)
% data is (T, nvar)

[T,n] = size(data);

% Get VC matrix
sigma = 1/T*data'*data;

[V,lamb] = eig(sigma); % diag(lamb) are the eigenvalues, V the eigenvectors

pc = data*V./n;
end
