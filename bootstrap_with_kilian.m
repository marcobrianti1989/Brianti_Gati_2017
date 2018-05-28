function [beta_tilde, data_boot2, beta_tilde_star,nonstationarities] ...
        = bootstrap_with_kilian(beta_hat, nburn, res, nsimul, ...
        which_correction, blocksize)

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
      [beta_hat_star1(:,:,i_simul),~,~] = reduform_var(data_boot1(:,:, i_simul), nlags);
end

psi_hat1 = mean(beta_hat_star1,3) - beta_hat;

% Step 1b
disp('Step 1b')
flag = test_stationarity(beta_hat');
if  flag == 0
      beta_tilde = beta_hat - psi_hat1;
else % nonstationary
      delt = 1;
      while flag == 1
            delt = delt - 0.01;
            beta_tilde = beta_hat - delt*psi_hat1;
            flag = test_stationarity(beta_tilde');
      end
end

%% 2nd bootstrap loop

% Step 2a
disp('Step 2a')
% Generate bootstrapped data samples
data_boot2 = data_boot(beta_tilde, nburn, res, nsimul, which_correction, blocksize);

beta_hat_star2  = zeros(nvar*nlags+1, nvar, nsimul);
beta_tilde_star = zeros(nvar*nlags+1, nvar, nsimul);

for i_simul = 1:nsimul
      [beta_hat_star2(:,:,i_simul),~,~] = reduform_var(data_boot2(:,:,i_simul), nlags);
end

% Step 2b
disp('Step 2b')
warning off
nonstationarities =zeros(nsimul,1);
for i_simul=1:nsimul
      flag = test_stationarity(beta_hat_star2(:,:,i_simul)');
      if  flag == 0
            beta_tilde_star(:,:,i_simul) = beta_hat_star2(:,:,i_simul) - psi_hat1; % we reuse psi_hat1 to save computational power
      else % nonstationary
            delt = 1;
            nonstationarities(i_simul,1) = 1;
            while flag == 1
                  if nonstationarities(i_simul,1) < 1000
                        delt = delt - 0.01;
                        beta_tilde_star(:,:,i_simul) = beta_hat_star2(:,:,i_simul) - delt*psi_hat1;
                        flag = test_stationarity(beta_tilde_star(:,:,i_simul)');
                        nonstationarities(i_simul,1) = nonstationarities(i_simul,1) + 1;
                  else
                        beta_tilde_star(:,:,i_simul) = beta_tilde;
                        flag = 0;
                        disp(['For simulation ' num2str(i_simul) ' we werent be able to reach stationarity.'])
                  end
            end
      end
end
      warning on