%Run a SVAR and identify news shocks and R&D shocks using Cholesky

% Marco Brianti, Laura Gáti, Oct 7 2017

clear all
close all

%Data Reading and Transformation
data = xlsread('dataset_23_sept_2017','Sheet1','B126:I283');
% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
TFP  = cumsum(data(:,1)); % TFP (levels in log)
RD   = log(data(:,2));  % R&D % this series was levels to start out with so don't cumsum <-- taking logs here induces stationarity of VAR - DISCUSS! If VAR nonstat and not cointegrated, estimation not possible.
IT   = log(data(:,5)); % IT investment; similar to RD
Mich = data(:,4); % the Mich index
GDP = log(data(:,6)); % real GDP % whether this guy's in logs or not doesn't seem to make a diff
C   = log(data(:,7)); % real cons % whether this guy's in logs or not doesn't seem to make a diff
H   = data(:,8); %hours worked

% Ordering in VAR
data_levels(:,1) = TFP;
data_levels(:,2) = H;
data_levels(:,3) = Mich;
data_levels(:,4) = IT;
data_levels(:,5) = GDP;
data_levels(:,6) = C;

% Generate automatically cell matrix of variable names for figures as well
% as define automatically which shocks to impose
which_shock = zeros(1,2);
varnames = cell(size(data_levels,2),1);
names = cell(1,2);
for i = 1:size(data_levels,2)
    if data_levels(:,i) == TFP
        varnames{i} = 'TFP';
    elseif data_levels(:,i) == H
        varnames{i} = 'H';
    elseif data_levels(:,i) == Mich
        varnames{i} = 'Mich index';
        which_shock(1,1) = i;
        names{1} = 'News shock';
    elseif data_levels(:,i) == IT
        varnames{i} = 'IT investment';
        which_shock(1,2) = i;
        names{2} = 'IT shock';
    elseif data_levels(:,i) == RD
        varnames{i} = 'R&D';
        which_shock(1,2) = i;
        names{2} = 'R&D shock';
    elseif data_levels(:,i) == GDP
        varnames{i} = 'GDP';
    elseif data_levels(:,i) == C
        varnames{i} = 'C';
    end
end

%Technical Parameters
max_lags   = 10;
nburn      = 0; %with the Kilian correction better not burning!!!
nsimul     = 5000; %5000
nvar       = size(data_levels,2);

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data_levels,max_lags);

%Run VAR imposing Cholesky
nlags = AIC;
[A,B,res,~] = sr_var(data_levels, nlags);

%Checking if the VAR is stationary
test_stationarity(B');

% Generate bootstrapped data samples
which_correction = 'blocks'; % [none, blocks] --> Choose whether to draws residuals in blocks or not.
q = 5;
dataset_boot = data_boot(B, nburn, res, nsimul, which_correction, q);

% Redo VAR nsimul times on the bootstrapped datasets
A_boot = zeros(nvar,nvar,nsimul);
B_boot = zeros(nvar*nlags+1,nvar,nsimul);
for i_simul = 1:nsimul
    [A_boot(:,:,i_simul), B_boot(:,:,i_simul), ~, ~] = ...
          sr_var(dataset_boot(:,:,i_simul), nlags);
end

% Kilian correction 
[B_corrected,  bias] = kilian_corretion(B, B_boot);
dataset_boot_corrected = data_boot(B_corrected, nburn, res, nsimul, which_correction, q);
A_boot_corrected = zeros(nvar,nvar,nsimul);
B_boot_corrected = zeros(nvar*nlags+1,nvar,nsimul);
for i_simul = 1:nsimul
    [A_boot_corrected(:,:,i_simul), B_boot_corrected(:,:,i_simul), ~, ~] = ...
          sr_var(dataset_boot_corrected(:,:,i_simul), nlags);
end
B_boot_test = mean(B_boot_corrected,3); %It should be very close to B
bias_test = sum(sum(abs(B - B_boot_test)));

A_boot = A_boot_corrected;
B_boot = B_boot_corrected;

%Calculate IRFs, bootstrapped CI and plot them
h=40; % horizon for IRF plots
sig = 0.90; % significance level
H = 100; % horizon for generation of IRFs
[IRFs, ub, lb] = genIRFs(A,A_boot,B,B_boot,H, sig);

plotIRFs(IRFs,ub,lb,h,which_shock, names, varnames)

% Variance decomposition
m = 40; %Horizon of the variance decomposition explained by the shocks
[vardec] = gen_vardecomp(IRFs,m,H);
[vardec_table] = vardecomp_table(vardec,which_shock,varnames,names);

%------------------------------------------------------------------------------------------


