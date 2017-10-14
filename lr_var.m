function [A,B,res,sigma] = lr_var(dataset, nlags)

% Rearrange data to run a VAR, outputting A, the impact matrix identified using LR restriction a la Blanchard & Quah.
% res = residuals from VAR estimation

[B,res,sigma] = rf_var(dataset, nlags);
nvar = size(B,2);

if nvar ==1
    A= 0;
else
    %LR restriction: A = C(1)^(-1) * [chol(C(1)*sigma* C(1)')]'  where
    %C(1) = [I - B]*(-1), the LR effect
    % Let G = [I - B]
    G = eye(nvar); % Need to build G with a loop because it is (nvar x nvar) while beta is ((nlags*nvar+1) x nvar)
    % Remove constant:
    B = B(2:end,:);
    
    B = B'; % I do this because I will fill up G always subtractin blocks of B corresponding to each lag
    for i = 1:nlags
        G = G - B(:, (i-1)*nvar+1:i*nvar);
    end
    C1 = inv(G);
    A  = inv(C1) * (chol(C1*sigma*C1'))';
    
end