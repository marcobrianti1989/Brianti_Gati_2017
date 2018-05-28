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
shocknames = {'News Shock','IT Shock'};
varnames

%Technical Parameters
max_lags        = 10;
nvar            = size(data,2);
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 200; %5000
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

return

% Bootstrap
which_correction = 'blocks'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks
[beta_tilde, data_boot2, beta_tilde_star, nonstationarities] = ...
      bootstrap_with_kilian_VECM(A', nburn, res', ...
      nsimul, which_correction, blocksize,r);

for i_simul = 1:nsimul
      %Run reduced form VECM
      %[alph_hat2, bet_hat2,Gam_hat2,sigma2] = RRR(data, nlags, r);
      [alph_hat_boot(:,:,i_simul),bet_hat_boot(:,:,i_simul),Pi_boot(:,:,i_simul),...
            Gam_hat_boot(:,:,i_simul),res_boot(:,:,i_simul),sigma_boot(:,:,i_simul)] = ...
            redu_VECM(data_boot2(:,:,i_simul), nlags, r); %It seems to be correct...
      %compared to a different method.
end

%For simplicity I keep the same impact. 
fake_impact = -[zeros(nvar,nvar-1) impact];
fake_impact_boot = zeros(nvar,nvar,nsimul);
for i_simul = 1:nsimul
      fake_impact_boot(:,which_shock,i_simul) = impact(:,:,i_simul);
end

return

return
close all
H = 40;
sig1 = 0;
sig2 = 0;
[IRFs, ub1, lb1, ub2, lb2] = genIRFs_VECM(fake_impact,0,A',0,H,sig1,sig2);
print_figs = 'no';
which_shock = [4];
simple_plotIRFs(IRFs,ub1,lb1,H,which_shock, print_figs)


%Creating and Printing figures
comment = [which_ID '_' char(varnames(5)) '_LR_hor_' num2str(LR_hor)];
print_figs = 'no';

[IRFs, ub1, lb1, ub2, lb2] = genIRFs(fake_impact,fake_impact_boot,...
      B,beta_tilde_star,H,sig1, sig2);

% % With Barsky & Sims-type ID, since you do a abs max, there are two
% % solutions: gam and -gam. So choose the one that makes sense. I'm doing
% % that so that news and IT have + effects on TFP.
if sum(IRFs(1,:,pos_IT)) < 0 % if the majority of TFP response is negative
      IRFs(:,:,pos_IT)   = -IRFs(:,:,pos_IT);
      ub1(:,:,pos_IT)     = - ub1(:,:,pos_IT);
      lb1(:,:,pos_IT)     = - lb1(:,:,pos_IT);
      ub2(:,:,pos_IT)     = - ub2(:,:,pos_IT);
      lb2(:,:,pos_IT)     = - lb2(:,:,pos_IT);
end
if sum(IRFs(1,:,pos_news)) < 0 % if the majority of TFP response is negative
      IRFs(:,:,pos_news)   = -IRFs(:,:,pos_news);
      ub1(:,:,pos_news)     = - ub1(:,:,pos_news);
      lb1(:,:,pos_news)     = - lb1(:,:,pos_news);
      ub2(:,:,pos_news)     = - ub2(:,:,pos_news);
      lb2(:,:,pos_news)     = - lb2(:,:,pos_news);
end

%Printing/Showing IRFs
h = 40;
% plotIRFs(IRFs,ub,lb,40,which_shocks,shocknames,varnames,which_ID,print_figs)
% plot_single_IRFs(IRFs,ub1,lb1,h,which_shocks,shocknames, varnames, which_ID, print_figs)
use_current_time = 0; % don't save the time
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,h,which_shocks,shocknames, varnames, '_', print_figs, use_current_time)

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

%Saving in Tex format the Variance Decomposition Matrix
fev_matrix = {'News', 'IT', 'Total'};
fev_matrix(2,:) = {num2str(FEV_news), num2str(FEV_IT), num2str(FEV_opt)};
disp('% of FEV of TFP explained:')
fev_matrix

export_FEV_matrix = 'no';
if strcmp(export_FEV_matrix,'yes') ==1
      fev_matrix_out = [FEV_news, FEV_IT, FEV_opt];
      rowLabels = {'Share of TFP FEV explained'};
      columnLabels = {'News', 'IT', 'Total'};
      matrixname   = 'FEVs';
      invoke_matrix_outputting(fev_matrix_out,matrixname,rowLabels,...
            columnLabels,comment);
end

toc
disp(varnames)
disp(datestr(now))

% Counterfactual
[s, IR] = get_structural_shocks_Forni(A,gam_opt,res);

which_shock_zero = 3; % set negative realizations of the IT shock to 0
y0 = data(1,:); % initial values of data
counterfactuals = counterfactual(s,IR,B,y0,which_shock_zero);

time = datetime(1989,7,1) + calquarters(0:size(s,1)); %ALWAYS CHECK IF THE INITIAL PERIOD IS THE SAME

figure
set(gcf,'color','w'); % sets white background color
set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
plot(time, data(:,1), 'b', 'linewidth', 2)
hold on
plot(time, counterfactuals(:,1), 'r', 'linewidth', 2)
legend('Original TFP series', 'Counterfactual TFP','Location','Southeast')
title('Shutting off negative IT productivity shocks')
set(gca, 'FontSize', 28)
grid on

save_counterfactual_fig = 'yes';
switch save_counterfactual_fig
      case 'yes'
            invoke_export_fig('counterfactual','',0)
end

% % Historical decompositions of TFP from the 2 structral shocks
% % IT:
% hd_IT = historical_decomposition(s(:,pos_IT),fake_impact,B, pos_IT, 1);
% % News:
% hd_news = historical_decomposition(s(:,pos_news),fake_impact,B, pos_news, 1);
% % sum of the two:
% hd = hd_IT + hd_news;
%
% % Get counterfactual when shutting off IT shocks 2000-Q3 onward
% IT_bubble = find(time=='01-Jul-2000');
% s_alternative = s;
% s_alternative(IT_bubble:end,3) = 0;
% % Redo counterfactuals for alternative scenario
% % IT:
% hd_IT_alt = historical_decomposition(s_alternative(:,pos_IT),fake_impact,B, pos_IT, 1);
% % News:
% hd_news_alt = historical_decomposition(s_alternative(:,pos_news),fake_impact,B, pos_news, 1);
% % sum of the two:
% hd_alt = hd_IT_alt + hd_news_alt;

[structural_shock_RyanID, ~] = ...
      get_structural_shocks_general(A,gam_opt,res,which_shocks);

