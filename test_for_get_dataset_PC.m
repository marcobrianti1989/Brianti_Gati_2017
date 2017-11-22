% Load and prepare all the series needed for Principal Components Analysis


clear
close all

tic % This takes about 50 sec.
% get_dataset_test_pc; % this downloads all the FRED series (81 series) and save them to individual
% % excel files

% Load in the SP500 series and convert them to quarterly
% Note: the monthly series starts in Jan 1950, so the quarterly starts in
% 1950-Q1.
base_path = pwd;
if exist([base_path '\Data\Single_series'], 'dir')
    cd([base_path '\Data\Single_series']) %for Microsoft
else
    cd([base_path '/Data/Single_series']) %for Mac
end

% [num, txt] = xlsread('YAHOO_SP500.csv','YAHOO_SP500.csv'); % xlsread
% can't read csv apparently... 
SP500 = csvread('YAHOO_SP500.csv',1,1);
cd(base_path)

closing_monthly = SP500(:,4);
closing_quarterly = closing_monthly(3:3:end,1); 
SP500 = closing_quarterly;
starting_date_SP500 = '1950-01-01';



toc

