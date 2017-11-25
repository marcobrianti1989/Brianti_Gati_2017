function [data, var_names] = read_data2(filename, sheet, range)
% This file is completely general function for reading in data. The only
% requirement for its use is that the dataset in Excel is structured in a
% specific way, which is:
% In one of the sheets, which you specify here using the sheet argument of
% this function, you should have the following rows:
% Row 1: Name of the variable (shorthand)
% Row 2: Transformation that is to be done on the variable (see possible transformations below)
% Row 3: In or Out: 1 if the series is to be used in the dataset (in a VAR
% maybe), 0 otherwise
% Row 4: No. series = just an accounting of how many series we have
% Row 5: Position: the position of the series in the dataset to be used (could play a role for Cholesky VARs for ex.)


% Create the correct path 
base_path = pwd;
if exist([base_path '\Data'], 'dir')
    addpath([base_path '\Data']) %for Microsoft
else
    addpath([base_path '/Data']) %for Mac
end

% Reading Excel
[numbers, variables]   = xlsread(filename,sheet,range);
variables              = variables(1,:); %Strings with the names of the variables
transformations        = numbers(1,:); %Numbers which represents which transformation to apply
In_or_out              = numbers(2,:); %Zero: variable out of the system, one: variables inside
nseries                = numbers(3,:); %Just the numner of the series (probably not useful)
position               = numbers(4,:); % the position of the variable in the VAR
full_data              = numbers(6:end,:); %Start from 6 since we avoid a useless initial NaN

nvar = size(variables,2);
% Transformations
% if 1: cumulate the series
% if 2: take the level
% if 3: take the log
% if 4: take the difference
% if 5: take the growth rate (diff after logs)
% else: keep the level but warning
for i_var = 1:nvar
      if     transformations(i_var) == 1
            full_data(:,i_var) = cumsum(full_data(:,i_var));
      elseif transformations(i_var) == 2 
            full_data(:,i_var) = full_data(:,i_var);
      elseif transformations(i_var) == 3
            full_data(:,i_var) = log(full_data(:,i_var));
      elseif transformations(i_var) == 4 
            full_data(:,i_var) = vertcat(nan, diff(full_data(:,i_var)));
      elseif transformations(:,i_var) == 5
            full_data(:,i_var) = vertcat(nan, diff(log(full_data(:,i_var))));
      else
            warning('Some variables may have a mispecified transformation')
      end     
end

% Selecting the VAR system we want to work with
data     = zeros(size(full_data,1), sum(In_or_out));
var_names= cell(1, sum(In_or_out));
for i_var = 1:nvar
      if In_or_out(i_var) == 1
            data(:,position(i_var)) = full_data(:,i_var);
            var_names(position(i_var)) = variables(i_var);
      elseif In_or_out(i_var) == 0
            % keep variable out of VAR
      else
            warning('The selection procedure of the system may be mispecified')
      end    
end

% Generalized Truncation
n_var_system = size(data,2);
threshold = -1/eps;
for i_var_system = 1:n_var_system
      loc(i_var_system) = find(data(:,i_var_system) > threshold, 1);
end
truncation_point = max(loc);
data = data(truncation_point:end,:);


end