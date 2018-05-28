function [beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
      = bootstrap_with_kilian_VECM(beta_hat, nburn, res, nsimul, ...
      which_correction, blocksize,r)

% This function does the entire bootstrap procedure using Kilian's
% correction, see Kilian 1998, p. 220.
% Inputs:
% beta_hat = the point estimate of beta from reduform_var
% Outputs:
% beta_tilde        = the corrected point estimate (nvar*nlags +1, nvar)
% data_boot2        = the generated simulated sample from the 2nd bootstrap loop (T-nlags,nvar,nsimul)
% beta_tilde_star   = bootstrapped betas (nvar*nlags +1, nvar, nsimul)
% nonstationarities = a vector that shows which simulations result in
% nonstationary betas and also tells you how many times they had to be
% discounted by delta until they became stationary.

%% 1st bootstrap loop

% Step 1a
disp('Step 1a')
% Generate bootstrapped data samples
data_boot1 = data_boot(beta_hat, nburn, res, nsimul, which_correction, blocksize);

[~, nvar, nsimul] = size(data_boot1);
if size(beta_hat,1)/size(beta_hat,2) == floor(size(beta_hat,1)/size(beta_hat,2))
      nlags = (size(beta_hat,1))/nvar;
else
      nlags = (size(beta_hat,1)-1)/nvar;
end

beta_hat_star1 = zeros(nvar*nlags+1, nvar, nsimul);
for i_simul = 1:nsimul
      [alph_hat_star1(:,:,i_simul),bet_hat_star1(:,:,i_simul),...
            Pi_star1(:,:,i_simul),Gam_hat_star1(:,:,i_simul),...
            U_star1(:,:,i_simul),sigma_star1(:,:,i_simul)] = ...
            redu_VECM(data_boot1(:,:, i_simul), nlags,r);
      %Setting the beta_hat (=A) matrix to obtain a SVAR functional form as follows
      % y(t) = A1*y(t-1) + ... + Ap*y(t-p) + B*eps(t)
      % Again, I am followinf Lutkepohl (2005)
      beta_hat_star1 = zeros(nvar,nvar*(nlags+1),nsimul);
      beta_hat_star1(:,1:nvar) = ...
            alph_hat_star1(:,:,i_simul)*bet_hat_star1(:,:,i_simul)' ...
            + eye(nvar) + Gam_hat_star1(:,1:nvar,i_simul);
      if nlags >= 2
            for i_lags = 1:nlags-1
                  beta_hat_star1(:,(i_lags*nvar)+1:(i_lags+1)*nvar) = ...
                        Gam_hat_star1(:,(i_lags*nvar)+1:(i_lags+1)*nvar,i_simul) ...
                        - Gam_hat_star1(:,((i_lags-1)*nvar)+1:i_lags*nvar,i_simul);
            end
      end
      beta_hat_star1(:,(nlags*nvar)+1:nvar*(nlags+1)) = ...
            - Gam_hat_star1(:,((nlags-1)*nvar)+1:nlags*nvar,i_simul);
end

psi_hat1 = mean(beta_hat_star1,3) - beta_hat;

% Step 1b
disp('Step 1b')
beta_tilde = beta_hat - psi_hat1;

%% 2nd bootstrap loop

% Step 2a
disp('Step 2a')
% Generate bootstrapped data samples
data_boot2 = data_boot(beta_tilde, nburn, res, nsimul, which_correction, blocksize);

beta_hat_star2  = zeros(nvar*nlags+1, nvar, nsimul);
beta_tilde_star = zeros(nvar*nlags+1, nvar, nsimul);

for i_simul = 1:nsimul
      [alph_hat_star2(:,:,i_simul),bet_hat_star2(:,:,i_simul),...
            Pi_star2(:,:,i_simul),Gam_hat_star2(:,:,i_simul),...
            U_star2(:,:,i_simul),sigma_star2(:,:,i_simul)] = ...
            redu_VECM(data_boot2(:,:,i_simul),nlags,r);
      %Setting the beta_hat (=A) matrix to obtain a SVAR functional form as follows
      % y(t) = A1*y(t-1) + ... + Ap*y(t-p) + B*eps(t)
      % Again, I am followinf Lutkepohl (2005)
      beta_hat_star2 = zeros(nvar,nvar*(nlags+1));
      beta_hat_star2(:,1:nvar) = alp*bet' + eye(nvar) + Gam(:,1:nvar);
      if nlags >= 2
            for i_lags = 1:nlags-1
                  beta_hat_star2(:,(i_lags*nvar)+1:(i_lags+1)*nvar) = ...
                        Gam(:,(i_lags*nvar)+1:(i_lags+1)*nvar) ...
                        - Gam(:,((i_lags-1)*nvar)+1:i_lags*nvar);
            end
      end
      beta_hat_star2(:,(nlags*nvar)+1:nvar*(nlags+1)) = - Gam(:,((nlags-1)*nvar)+1:nlags*nvar);
end

% Step 2b
disp('Step 2b')

for i_simul=1:nsimul
      beta_tilde_star(:,:,i_simul) = beta_hat_star2(:,:,i_simul) - psi_hat1;
end

