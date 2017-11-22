% Load and prepare all the series needed for Principal Components Analysis


clear
close all

tic % This takes about 50 sec.
get_dataset_test_pc; % this downloads all the FRED series (81 series) and save them to individual
% excel files

%% Load in the SP500 series and convert them to quarterly
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
closing_monthly = SP500(:,4);
closing_quarterly = closing_monthly(3:3:end,1); 
SP500 = closing_quarterly;
name_SP = 'Quarterly_SP500';
starting_date_SP500 = '_Starts_1950_q1';
date_aggregated = datestr(today);
comment = ['Date aggregated: ', date_aggregated];
filename = [name_SP, starting_date_SP500];
xlswrite(filename,SP500)


%% Download Michigan index business conditions series as a CSV file in Single series

url_mich_bus_csv = 'https://data.sca.isr.umich.edu/sda-public/tmpdir/AAc51CRa.csv'; % business conditions
filename_mich = 'mich'; % this is the name of the CSV file
outfilename = websave(filename_mich,url_mich_bus_csv);

%% Download Fernald's TFP database

url_fernald = 'https://www.frbsf.org/economic-research/files/quarterly_tfp.xlsx';
web_tfp = webread(url_fernald);
filename_fernald = 'fernald';
outfilename = websave(filename_fernald,url_fernald);


cd(base_path)
toc

