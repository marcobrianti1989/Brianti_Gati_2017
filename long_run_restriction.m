function [B0] = long_run_restriction(beta, omega)

nvar = size(beta,2);
nlag = size(beta,1)/nvar;

beta = beta';

G = eye(nvar);
for i = 1:nlag
    G = G - beta(:, (i-1)*nvar+1:i*nvar);
end

C = inv(G);

B0 = inv(C)*(chol(C*omega*C')');