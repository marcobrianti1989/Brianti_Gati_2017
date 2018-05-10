function [beta_tilde, data_boot2, beta_tilde_star] ...
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
nlags = (size(beta_hat,1)-1)/nvar - 1;

%beta_hat_star1 = zeros(nvar*nlags+1, nvar, nsimul);
for i_simul = 1:nsimul
      [alph_hat_star1(:,:,i_simul),bet_hat_star1(:,:,i_simul),...
            Pi_star1(:,:,i_simul),Gam_hat_star1(:,:,i_simul),...
            res_star1(:,:,i_simul),sigma_star1(:,:,i_simul)] ...
            = redu_VECM(data_boot1(:,:, i_simul), nlags, r);
      %Run structural VECM just for beta
      % Notation note: here we call beta the estimated AR coefficient that in the
      % VECM literature a la Lutkepohl is usually denoted by A.
      [~, ~, beta_hat_star1(:,:,i_simul)] ...
            = structural_VECM(alph_hat_star1(:,:,i_simul),bet_hat_star1(:,:,i_simul),...
            Gam_hat_star1(:,:,i_simul),res_star1(:,:,i_simul),sigma_star1(:,:,i_simul),nlags,r);
end

beta_hat = beta_hat(2:end,:);
psi_hat1 = (mean(beta_hat_star1,3))' - beta_hat;
beta_tilde = beta_hat - psi_hat1;
%Here we are not controlling for stationarity since the VECM is non
%stationary by construction.

%% 2nd bootstrap loop

% Step 2a
disp('Step 2a')
% Generate bootstrapped data samples
beta_tilde = [zeros(size(beta_tilde,2),1)'; beta_tilde];
data_boot2 = data_boot(beta_tilde, nburn, res, nsimul, which_correction, blocksize);

%beta_hat_star2  = zeros(nvar*nlags+1, nvar, nsimul);
beta_tilde_star = zeros(nvar*nlags+1, nvar, nsimul);

for i_simul = 1:nsimul
     % [beta_hat_star2(:,:,i_simul),~,~] = reduform_var(data_boot2(:,:,i_simul), nlags);
            [alph_hat_star2(:,:,i_simul),bet_hat_star2(:,:,i_simul),...
            Pi_star2(:,:,i_simul),Gam_hat_star2(:,:,i_simul),...
            res_star2(:,:,i_simul),sigma_star2(:,:,i_simul)] ...
            = redu_VECM(data_boot2(:,:, i_simul), nlags, r);
            %Run structural VECM just for beta
      [~, ~, beta_hat_star2(:,:,i_simul)] ...
            = structural_VECM(alph_hat_star2(:,:,i_simul),bet_hat_star2(:,:,i_simul),...
            Gam_hat_star2(:,:,i_simul),res_star2(:,:,i_simul),sigma_star2(:,:,i_simul),nlags,r);
end

% Step 2b
disp('Step 2b')
beta_tilde_star = permute(beta_hat_star2,[2 1 3]) - psi_hat1;