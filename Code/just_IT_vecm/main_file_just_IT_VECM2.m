% This file corresponds to main_file_VECM, it is a tool to test if the
% function VECM is working

% Marco Brianti, Laura Gáti, Mar 13 2018
clear
close all

current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

if exist([base_path '\Data'], 'dir')
      addpath([base_path '\Data']) %for Microsoft
else
      addpath([base_path '/Data']) %for Mac
end

%Data Reading and Transformation
filename = 'dataset_main';
sheet    = 'Data';
range    = 'B1:W287';
[data, varnames] = read_data2(filename, sheet, range);
shocknames = {'IT Shock','Fake'};
varnames

%Technical Parameters
max_lags        = 10;
nvar            = size(data,2);
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 1000; %5000
sig1            = 0.95; % significance level
sig2            = 0.9; % a 2nd sig. level
H               = 40; %40; % horizon for generation of IRFs
h               = H; %40; % horizon for IRF plots
which_variable  = 1; % select TFP as the variable whose FEV we wanna max

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data,max_lags);
if AIC >= 4
      nlags = 1; % 2;
      warning(['AIC > 4, setting nlags to ', num2str(nlags) ])
else
      nlags = AIC;
end

%Run reduced form VECM
r = 1; %number of cointegrating vectors
%[alph_hat2, bet_hat2,Gam_hat2,sigma2] = RRR(data, nlags, r);
[alph_hat,bet_hat,Pi,Gam_hat,res,sigma] = redu_VECM(data, nlags, r); %It seems to be correct...
%compared to a different method.

%Run structural VECM
diff_BBp_sigma = 1;
j = 0;
while diff_BBp_sigma >= 10^(-6)
      j = j + 1
      [B, Xi, A, diff_BBp_sigma] = structural_VECM(alph_hat,bet_hat,Gam_hat,res,sigma,nlags,r);
      diff_BBp_sigma
end

% Implement the "just IT" ID strategy in a VAR
which_variable = 2;
which_shock = 4;
[impact, impact_IT_opt, gam_opt]  = ...
      just_IT_ID_VECM(which_variable,which_shock,B,Xi,r);

% Bootstrap
which_correction = 'blocks'; %[none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks
% [beta_tilde, data_boot2, beta_tilde_star, nonstationarities] = ...
%       bootstrap_with_kilian_VECM(A', nburn, res', ...
%       nsimul, which_correction, blocksize,r);

data_boot2 = ...
      data_boot(A', nburn, res', nsimul, which_correction, blocksize);

A_boot = zeros(nvar,nvar*(nlags+1),nsimul);
for i_simul = 1:nsimul
      %Run reduced form VECM
      %[alph_hat2, bet_hat2,Gam_hat2,sigma2] = RRR(data, nlags, r);
      [alph_hat_boot(:,:,i_simul),bet_hat_boot(:,:,i_simul),Pi_boot(:,:,i_simul),...
            Gam_hat_boot(:,:,i_simul),res_boot(:,:,i_simul),sigma_boot(:,:,i_simul)] = ...
            redu_VECM(data_boot2(:,:,i_simul), nlags, r); %It seems to be correct...
      %compared to a different method.
      A_boot(:,1:nvar,i_simul) = ...
            alph_hat_boot(:,:,i_simul)*bet_hat_boot(:,:,i_simul)' ...
            + eye(nvar) + Gam_hat_boot(:,1:nvar,i_simul);
      if nlags >= 2
            for i_lags = 1:nlags-1
                  A_boot(:,(i_lags*nvar)+1:(i_lags+1)*nvar,i_simul) = ...
                        Gam_hat_boot(:,(i_lags*nvar)+1:(i_lags+1)*nvar,i_simul) ...
                        - Gam_hat_boot(:,((i_lags-1)*nvar)+1:i_lags*nvar,i_simul);
            end
      end
      A_boot(:,(nlags*nvar)+1:nvar*(nlags+1),i_simul) = ...
            - Gam_hat_boot(:,((nlags-1)*nvar)+1:nlags*nvar,i_simul);
end
A_boot = permute(A_boot,[2 1 3]);

%For simplicity I keep the same impact.
fake_impact = -[zeros(nvar,nvar-1) impact];
fake_impact_boot = zeros(nvar,nvar,nsimul);
for i_simul = 1:nsimul
      fake_impact_boot(:,end,i_simul) = -impact;
end

H = 40;
sig1 = .90;
sig2 = .95;
[IRFs, ub1, lb1, ub2, lb2] = ...
      genIRFs_VECM(fake_impact,fake_impact_boot,A',A_boot,H,sig1,sig2);
print_figs = 'no';
which_shock = [4];

%Creating and Printing figures
print_figs = 'no';
base_path = 'C:\Users\th3\Documents\BG_2017';
which_ID = 'Just_IT_VECM';
use_current_time = 0; % don't save the time
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,H,...
      which_shock,shocknames, varnames,which_ID, print_figs, ...
      use_current_time,base_path)

return

%Forni&Gambetti Orthogonality Test
do_FG_test = 'no';
switch do_FG_test
      case 'yes'
            filename_PC       = 'Dataset_test_PC';
            sheet_PC          = 'Quarterly';
            range_PC          = 'B2:DC287';
            first_n_PCs       = 10;
            mlags             = 1;
            [pvalue_news_shock, pvalue_IT_shock] = ...
                  Forni_Gambetti_orthogonality_test(filename_PC,...
                  sheet_PC,range_PC,first_n_PCs,A,gam_opt,res,which_shocks,mlags)
end

% Build the denominator of variance decomposition
N = null(gam_opt');
D_null = [N gam_opt];
D_null*D_null'; % this is just a check, should be the identity.
impact_vardec = B*D_null; % where A is the chol.
[IRFvardec, ~, ~, ~, ~] = ...
      genIRFs_VECM(impact_vardec,0,A',0,H,sig1,sig2);

vardec = gen_vardecomp(IRFvardec,h,h);
shareIT_on_TFP = vardec(1,end); % "end" b/c we put gam_opt as last.

%Get Structural Shocks
[s_shock_just_IT, ~] = ...
      get_structural_shocks_general(B,gam_opt,res',which_shock);




