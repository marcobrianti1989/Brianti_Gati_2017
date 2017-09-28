%% Run a SVAR and identify news shocks and R&D shocks

% Marco Brianti, Laura Gáti, Sep 28 2017

clear all
close all

data = xlsread('dataset_23_sept_2017','Sheet1','B126:F283');
% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)A
data_levels(:,1) = cumsum(data(:,1));
data_levels(:,2) = data(:,5); % this series was levels to start out with so we take logs
% 5 means we take IT investment instead of R&D
% data_levels(:,3) = cumsum(data(:,3)); %ignore this since it's SPF
data_levels(:,3) = data(:,4); % the Mich index take the third column in data_levels.

%Wrong order in the dataset! Rearrange to have: RD is the last! TFP first and SPF second!
y      = zeros(size(data_levels,1),size(data_levels,2));
y(:,1) = data_levels(:,1);
y(:,2) = data_levels(:,3);
y(:,3) = data_levels(:,2);

data_levels = y;

lag_number = 2;
which_shock = 2;
total_extractions = 5000;

tic
[A, B, B_boot, B_boot_Kilian, mshock, A99, A1, A95, A5, A86, A16] = ...
    cholboot(data_levels,lag_number,which_shock,total_extractions);
toc

% For checking purposes: an average beta_boot
avg_B_boot = mean(B_boot_Kilian,3);

%% Redo stuff with our new shiny simple (and quick) code
nlags = lag_number;
nburn = 200;
nsimul = 500;
nvar = size(data_levels,2);

% Run VAR imposing Cholesky
[A2,B2,res] = sr_var(data_levels, nlags);

% Generate bootstrapped data samples
dataset_boot = data_boot(B2, nburn, res, nsimul); % <--- TO DO: draw shocks in blocks

% Redo VAR nsimul times on the bootstrapped datasets
A_boot = zeros(nvar,nvar,nsimul);
B_boot = zeros(nvar*nlags+1,nvar,nsimul);
for i=1:nsimul
    [A_boot(:,:,i), B_boot(:,:,i), ~] = sr_var(dataset_boot(:,:,i), nlags);
end

average_B_boot = mean(B_boot,3);
average_A_boot = mean(A_boot,3);
% TO DO: Kilian correction


% Calculate IRFs, bootstrapped CI and plot them
h=80;
which_shock = [2 3];
names = {'News shock', 'R&D shock'}; % shock names in order of appearance
varnames = {'TFP','Mich index', 'R&D'} ; % variable names in order of appearance
sig = 0.9; % significance level
plotIRFs(A2,A_boot,B,B_boot,h,which_shock, names, varnames,sig);