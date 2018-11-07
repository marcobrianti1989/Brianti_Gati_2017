% main_file_just_IT.m  - implements a SVAR to identify the effect of an IT
% shock alone as a shock that maximizes FEV IT on impact and has 0 impact
% effect on TFP
%
% Code by Marco Brianti and Laura Gati, Boston College, 2018
%**************************************************************

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

tic

%Data Reading and Transformation
filename = 'Quarterly';
sheet    = 'Quarterly Data';
range    = 'B1:V275';
warning off
[data, varnames] = read_data2(filename, sheet, range);
warning on
shocknames = {'ICT Shock'};

for i = 1:size(data,2)
      eval([varnames{i} ' = data(:,i);'])
end

% Proper Transformations - All the variables should be in logs
percapita = 1;
if percapita == 1
      Hours               = Hours + Employment - Population; %Average weekly hours over population
      Consumption         = NonDurableCons + ServiceCons - Population;
      Investment          = Investment + DurableCons - Population;
      GDP                 = GDP - Population;
      SP5001              = SP5001 - Population - GDPDef;
      SP5002              = SP5002 - Population - GDPDef;
      ICTInvestment       = RealICTInvestment - Population;
else
      Consumption         = NonDurableCons + ServiceCons;
      Investment          = Investment + DurableCons;
      SP5001              = SP5001 - GDPDef;
      SP5002              = SP5002 - GDPDef;
      ICTInvestment       = RealICTInvestment;
end
which_shock = 2; %[pos_IT];

%Technical Parameters
max_lags        = 10;
nburn           = 0; %with the Kilian correction better not burning!!!
nsimul          = 50; %5000
nvar            = size(data,2);
sig1            = 0.9; % significance level
sig2            = 0.95; % a 2nd sig. level
H               = 40; %40; % horizon for generation of IRFs
h               = 40; %40; % horizon for IRF plots
which_variable  = which_shock; % select IT as the variable whose FEV we wanna max

% Define the system1
system_names  = {'TFPUtil','ICTInvestment','GDP','Consumption',...
      'Investment','Hours','RelativePriceICT'};
for i = 1:length(system_names)
      system(:,i) = eval(system_names{i});
end
which_variable = find(strcmp('ICTInvestment', system_names));

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
nlags = 4;

%Run VAR imposing Cholesky
[A,B,res,sigma] = sr_var(system, nlags);
%Checking if the VAR is stationary
test_stationarity(B');
% Get Chol Structural Shocks
ss          = (inv(A)*res')';
corr(ss);
ssICT       = ss(:,2);
ssICT_time  = [Time(1+nlags:end) ssICT]; 
save('ssICT_time','ssICT_time')
plot(Time(1+nlags:end),ssICT)
asdf
% Implement the "just IT" ID strategy in a VAR
H_max = 40;
%[impact, impact_IT_opt, gam_opt]  = just_IT_ID(which_variable,which_shock,A);

% Bootstrap
which_ID = 'just_IT';
which_correction = 'none'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
blocksize = 5; % size of block for drawing in blocks
[beta_tilde, data_boot2, beta_tilde_star, nonstationarities] = ...
      bootstrap_with_kilian(B, nburn, res, ...
      nsimul, which_correction, blocksize);


% Get "bootstrapped A" nsimul times
B_boot = zeros(size(B,1),size(B,2),nsimul);
A_boot = zeros(size(A,1),size(A,2),nsimul);
for i_simul=1:nsimul
      [A_boot(:,:,i_simul), B_boot(:,:,i_simul),~,~] = sr_var(data_boot2(:,:,i_simul), nlags);
      % Get bootstrapped confidence intervals nsimul times
      disp(['Iteration ' num2str(i_simul) ' out of ' num2str(nsimul)])
      %[impact_boot(:,i_simul), ~, ~]  = just_IT_ID(which_variable,which_shock,A_boot);
end

%Creating a fake matrix for the IRF of the point estimation
% fake_impact = zeros(size(system,2),size(system,2));
% fake_impact(:,which_shock) = impact;
% % Create a fake matrix for the IRFs of the bootstrap
% fake_impact_boot = zeros(size(system,2),size(system,2),nsimul);
% for i_simul = 1:nsimul
%       fake_impact_boot(:,which_shock,i_simul) = impact_boot(:,i_simul);
% end

%Creating and Printing figures
comment = [which_ID '_'];
print_figs = 'no';

[IRFs, ub1, lb1, ub2, lb2] = genIRFs(A,A_boot,...
      B,B_boot,H,sig1, sig2);

use_current_time = 0; % don't save the time
plot_single_IRFs_2CIs(IRFs,ub1,lb1,ub2,lb2,h,which_shock,shocknames,...
      system_names, 'empirical_noH', print_figs, use_current_time, base_path)

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
                  sheet_PC,range_PC,first_n_PCs,A,gam_opt,res,which_shock,mlags)
end

% Build the denominator of variance decomposition
% N = null(gam_opt');
% D_null = [N gam_opt];
% D_null*D_null'; % this is just a check, should be the identity.
% impact_vardec = A*D_null; % where A is the chol.
% [IRF_vardec, ~, ~, ~, ~] = genIRFs(impact_vardec,0,B,0,H, sig1, sig2);

vardec          = gen_vardecomp(IRF,h,h);
shareIT_on_TFP  = vardec(1,end); % "end" b/c we put gam_opt as last.
h1     = 1;
h4     = 4;
h8     = 8;
h16    = 16;
h24    = 24;
h40    = 40;
h_vec  = [h1 h4 h8 h16 h24 h40];
for ih = 1:length(h_vec)
      vardec = gen_vardecomp(IRF_vardec,h_vec(ih),h);
      table_vardec(:,ih) = vardec(:,end);
end


%Get Structural Shocks
[s_shock_just_IT, ~] = ...
      get_structural_shocks_general(A,gam_opt,res,which_shock);
%plot(s_shock_just_IT)

varnames

disp('stopping here')

return

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
