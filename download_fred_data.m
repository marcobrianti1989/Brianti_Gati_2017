function download_fred_data

%% Trying to read content from FRED specifying frequency - works
series_id = 'FEDFUNDS';
observation_start = '1947-01-01';
observation_end   = datestr(today,'yyyy-mm-dd');
units  = 'log';
frequency = 'q';
aggregation_method = 'eop';
[output] = getFredData(series_id, observation_start, observation_end, units, frequency, aggregation_method)
time_fedfunds = datestr(output.Data(:,1));