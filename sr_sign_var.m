function [A, B] = sr_sign_var(data_levels, nlags,sign_res)
% Implementation of Rubio, Waggoner, Zha 2009. Algorithm 2, p. 688
% Algorithm 2. Let (A0,A+) be any given value of the unrestricted structural parameters.
% (Step 1) Draw an independent standard normal n × n matrix X and let X = Q R be the QR
% decomposition of X with the diagonal of R normalized to be positive.
% (Step 2) Let P = Q and generate impulse responses from (A0*P)^(?1) and B = Aplus*A0^(?1)
% (Step 3) If these impulse responses do not satisfy the sign restrictions, return to Step 1.

% Output: A0*P which is the new impact matrix

sign_res = sign_res(:);
ind = find(sign_res~=0);
sign_res = sign_res(ind);

[T, nvar] = size(data_levels);
% Estimate a VAR where you obtain beta and an initial impact matrix A and
% almost-VC matrix sigma
[A,B,~,~] = sr_var(data_levels, nlags);

maxnum = 1000;
i = 1;
check = 1;
while i <= maxnum && check > 0
      % (1) Draw n x n random standard normal matrix
      X = normrnd(0,1,nvar,nvar);
      [Q,R]=qr(X); % Matlab function implementing the QR decomposition
      % Normalize (dubious in paper...)
      for j=1:nvar
            if R(j,j)<0
                  Q(:,j)=-Q(:,j);
                  R(j,:)=-R(j,:); % we don't actually need this, but we keep it here to be clear
            end
      end
      
      % (2) Define P and generate IRFs
      P = Q;
      A0P = (A*P)^(-1); % A0P^-1 is the new structural matrix candidate
      
      % %Calculate IRFs
      % sig = 0.90; % significance level
      % H = 100; % horizon for generation of IRFs
      % [IRFs, ~, ~] = genIRFs(A0P^(-1),0,B,0,H, sig);
      
      % (3) Check whether sign restrictions are satisfied (for now only on impact)
      A0P = A0P(:);
      actual_signs = sign(A0P);
      
      check = sum(abs(sign_res - actual_signs(ind)));
      i = i+1
end
A = reshape(A0P, nvar, nvar);




end

