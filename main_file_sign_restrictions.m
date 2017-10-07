%Run a SVAR and identify news shocks and R&D shocks using sign restrictions

% Marco Brianti, Laura Gáti, Oct 7 2017

clear all
close all

%Data Reading and Transformation
data = xlsread('dataset_23_sept_2017','Sheet1','B126:J283');
% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
TFP  = cumsum(data(:,1)); % TFP (levels in log)
RD   = log(data(:,2));  % R&D % this series was levels to start out with so don't cumsum <-- taking logs here induces stationarity of VAR - DISCUSS! If VAR nonstat and not cointegrated, estimation not possible.
IT   = log(data(:,5)); % IT investment; similar to RD
Mich = data(:,4); % the Mich index
GDP = log(data(:,6)); % real GDP % whether this guy's in logs or not doesn't seem to make a diff
C   = log(data(:,7)); % real cons % whether this guy's in logs or not doesn't seem to make a diff
H   = data(:,8); %hours worked
P   = vertcat(nan, diff(log(data(:,9)))); %price of IT goods
P(isnan(P)) = -999;

% Ordering in VAR
data_levels(:,1) = TFP;
data_levels(:,2) = H;
data_levels(:,3) = Mich;
data_levels(:,4) = IT;
data_levels(:,5) = GDP;
data_levels(:,6) = C;
data_levels(:,7) = P;

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
      elseif data_levels(:,i) == P
            varnames{i} = 'Price IT';
      end
end

% Truncate dataset to the shortest variable
P(P==-999) = nan;
start = find(isnan(P) < 1,1,'first');
data_levels = data_levels(start:end,:);

%Technical Parameters
max_lags   = 8;
nburn      = 0; %with the Kilian correction better not burning!!!
nsimul     = 5000; %5000
nvar       = size(data_levels,2);

%%Checking the number of lags over BIC, AIC, and HQ (see 'Lecture2M' in our folder)
[AIC,BIC,HQ] = aic_bic_hq(data_levels,max_lags);
nlags = AIC;

% Build sign restriction matrix
sign_res = zeros(nvar,nvar);
sign_res(end,3) = 1; % posiive sign restriction 
sign_res(end,4) = -1; % negative sign restriction

% Now impose sign restrictions
[A, B] = sr_sign_var(data_levels, nlags,sign_res);

%Calculate IRFs, bootstrapped CI and plot them
h=40; % horizon for IRF plots
sig = 0.90; % significance level
H = 100; % horizon for generation of IRFs
[IRFs, ub, lb] = genIRFs(A,0,B,0,H, sig);

plotIRFs(IRFs,ub,lb,h,which_shock, names, varnames)

