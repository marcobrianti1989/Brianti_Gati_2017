clear
close all

% Trying to read content from FRED specifying frequency - works
series_id = 'FEDFUNDS';
observation_start = '1947-01-01';
observation_end   = datestr(today,'yyyy-mm-dd');
units  = 'log';
frequency = 'q';
aggregation_method = 'eop';

dataset_PC = get_dataset_test_pc(series_id, observation_start, ...
      observation_end, units, frequency, aggregation_method, ondate, realtime_end);