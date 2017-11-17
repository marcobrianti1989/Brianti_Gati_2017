function [data_levels, shocknames, varnames, which_shock, q] = read_data(filename, sheet, range)
% A file for reading in data specifically for this project. TO DO: Make
% this general and clean up data reading in general.
base_path = pwd;
if exist([base_path '\Data'], 'dir')
    addpath([base_path '\Data']) %for Microsoft
else
    addpath([base_path '/Data']) %for Mac
end
data = xlsread(filename,sheet,range);
% Cumulate growth variables to levels (log levels to be precise, b/c growth
% rates are calculated as log diffs)
TFP   = cumsum(data(:,1)); % TFP (levels in log)
RD    = log(data(:,2));  % R&D % this series was levels to start out with so don't cumsum <-- taking logs here induces stationarity of VAR - DISCUSS! If VAR nonstat and not cointegrated, estimation not possible.
IT    = log(data(:,5)); % IT investment; similar to RD
Mich  = data(:,4); % the Mich index
GDP   = log(data(:,6)); % real GDP % whether this guy's in logs or not doesn't seem to make a diff
C     = log(data(:,7)); % real cons % whether this guy's in logs or not doesn't seem to make a diff
H     = data(:,8); %hours worked
P_IT  = log(data(:,9)); %price of IT goods 
P_IT  = vertcat(nan, diff(log(data(:,9)))); %price of IT goods (to be correct, this is inflation in price index of IT gods)
Pi    = log(data(:,10)); % log CPI
Pi = vertcat(nan, diff(log(data(:,10)))); % CPI inflation
rel_price = P_IT - Pi; % this is the ratio of IT price inflation over CPI inflation
rel_price(isnan(rel_price)) = -999;
P_K   = log(data(:,11)); % log price of capital
% P_K   = vertcat(-999, diff(P_K)); % try price of capital in differences

[rho, instrIT, ~] = quick_ols(IT(1:end-1,:), IT(2:end,:));
instrIT = vertcat(nan, instrIT);
% instrIT = vertcat(nan, IT(1:end-1,:));
instrIT(isnan(instrIT)) = -999;


% Ordering in VAR
data_levels(:,1) = TFP;
% data_levels(:,2) = H;
data_levels(:,2) = Mich;
data_levels(:,3) = IT; %RD;
% data_levels(:,4) = instrIT;
data_levels(:,4) = GDP;
data_levels(:,5) = C;
% data_levels(:,7) = RD;
data_levels(:,6) = rel_price;
% data_levels(:,6) = P_K;


% Generate automatically cell matrix of variable names for figures as well
% as define automatically which shocks to impose
which_shock = zeros(1,2);
varnames = cell(size(data_levels,2),1);
shocknames = cell(1,2);
do_truncation  = 'no';
do_truncation2 = 'no';
do_truncation3 = 'no';
do_truncation4 = 'no';
do_truncation5 = 'no';
q = NaN;
for i = 1:size(data_levels,2)
    if data_levels(:,i) == TFP
        varnames{i} = 'TFP';
    elseif data_levels(:,i) == H
        varnames{i} = 'H';
    elseif data_levels(:,i) == Mich
        varnames{i} = 'Mich index';
        which_shock(1,1) = i;
        shocknames{1} = 'News shock';
        do_truncation5 = 'yes';
    elseif data_levels(:,i) == IT
        varnames{i} = 'IT investment';
        which_shock(1,2) = i;
        shocknames{2} = 'IT shock';
    elseif data_levels(:,i) == RD
        varnames{i} = 'R&D';
        which_shock(1,2) = i;
        shocknames{2} = 'R&D shock';
    elseif data_levels(:,i) == GDP
        varnames{i} = 'GDP';
    elseif data_levels(:,i) == C
        varnames{i} = 'C';
    elseif data_levels(:,i) == P_IT
        varnames{i} = 'Price IT';
        do_truncation = 'yes';
    elseif data_levels(:,i) == Pi
        varnames{i} = 'CPI Inflation';
        do_truncation2 = 'yes';
    elseif data_levels(:,i) == instrIT
        varnames{i} = 'IT investment (IV)';
        shocknames{2} = 'IT shock';
        which_shock(1,2) = i;
        do_truncation3 = 'yes';
    elseif data_levels(:,i) == rel_price
        varnames{i} = 'Relative price of IT';
        do_truncation4 = 'yes';
        q = i; % position of rel. price of IT
    elseif data_levels(:,i) == P_K
        varnames{i} = 'Capital price';
        q = i; % pretend like this was the position of rel. price of IT (we're using P_K as a proxy for P_IT)
    end
end

if strcmp(do_truncation, 'yes')
    % Truncate dataset to the shortest variable
    P_IT(P_IT==-999) = nan;
    start = find(isnan(P_IT) < 1,1,'first');
    data_levels = data_levels(start:end,:);
    
elseif strcmp(do_truncation2, 'yes')
    % Truncate dataset to the shortest variable
    Pi(Pi==-99) = nan;
    start = find(isnan(Pi) < 1,1,'first');
    data_levels = data_levels(start:end,:);
    
elseif strcmp(do_truncation3, 'yes')
    % Truncate dataset to the shortest variable
    instrIT(instrIT==-999) = nan;
    start = find(isnan(instrIT) < 1,1,'first');
    data_levels = data_levels(start:end,:);
elseif strcmp(do_truncation4, 'yes')
    % Truncate dataset to the shortest variable
    rel_price(rel_price==-999) = nan;
    start = find(isnan(rel_price) < 1,1,'first');
    data_levels = data_levels(start:end,:);
elseif strcmp(do_truncation5, 'yes')
    % Truncate dataset to the shortest variable
    Mich(Mich==-999) = nan;
    start = find(isnan(Mich) < 1,1,'first');
    data_levels = data_levels(start:end,:);
end