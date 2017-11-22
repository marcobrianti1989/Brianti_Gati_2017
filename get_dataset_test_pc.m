function get_dataset_test_pc
% This function is automatically bulding a large dataset with aggregate US
% data which will then be used to performe the principal component
% analysis.

observation_start = '1947-01-01';
observation_end   = datestr(today,'yyyy-mm-dd');
units  = 'log';
frequency = 'q';
aggregation_method = 'eop';

series_id_vector = { 'GDPC1',	'GNPC96',	'NICUR',	'DPIC96',	'OUTNFB',	'FINSLC1',	'FPIC1',	'PRFIC1', ...
    'PNFIC1',	'GPDIC1',	'PCECC96',	'PCNDGC96',	'PCDGCC96',	'PCESVC96',	'GPSAVE',	'FGCEC1',	'FGEXPND', ...
    'FGRECPT', 'CBIC1',	'EXPGSC1',	'IMPGSC1',	'CP',	'NFCPATAX',	'CNCF',	'DIVIDEND',	'HOANBS',	'OPHNFB', ...
    'UNLPNBS',	'ULCNFB',	'WASCUR',	'COMPNFB',	'COMPRNFB',	'GDPCTPI',	'GNPCTPI',	'GDPDEF',	'GNPDEF', ...
    'INDPRO',	'IPBUSEQ',	'IPCONGD',	'IPDCONGD',	'IPFINAL',	'IPMAT',	'IPNCONGD',	'AWHMAN',	'AWOTMAN', ...
    'CIVPART',	'CLF16OV',	'CE16OV',	'USPRIV',	'USGOOD',	'SRVPRD',	'UNEMPLOY',	'UEMPMEAN',	'UNRATE', ...
    'HOUST',	'FEDFUNDS',	'TB3MS',	'GS1',	'GS10',	'AAA',	'BAA',	'MPRIME',	'M1SL',	'M2MSL',	'M2SL', ...
    'BUSLOANS',	'CONSUMER',	'LOANINV',	'REALLN',	'TOTALSL',	'CPIAUCSL',	'CPIULFSL',	'CPILEGSL',	'CPILFESL', ...
    'CPIENGSL',	'CPIUFDSL',	'PPICPE',	'PPICRM',	'PPIFCG',	'PPIFGS',	'OILPRICE'}; %

% data_fred  = zeros(283,size(series_id_vector,2));
for i = 1:size(series_id_vector,2)
    series_id = series_id_vector{i};
series = getFredData(series_id, observation_start, observation_end, units, frequency, aggregation_method);

base_path = pwd;
if exist([base_path '\Data\Matlab_download'], 'dir')
    cd([base_path '\Data\Matlab_download']) %for Microsoft
else
    cd([base_path '/Data/Matlab_download']) %for Mac
end



warning off
xlswrite(series_id,series.Data)

cd(base_path)
% length_this_series = size(series.Data,1);
% length_difference = 283-length_this_series;
% if series.Data(1,1) == datenum(observation_start) % the missing obs are at the end (i.e. the most recent ones)
%     series.Data(:,2) = 
% data_fred(:,i)  = series.Data(:,2);
warning on
end
% times_fred = series.Data(:,1); % Matlab datenums
% times_fred_string = datestr(times_fred,'yyyy-mm-dd');