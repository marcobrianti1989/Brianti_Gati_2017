function [data, var_names] = read_data2(filename, sheet, range)

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
full_data              = numbers(5:end,:); %Start from 5 since we avoid a useless initial NaN

nvar = size(variables,2);
% Transformations
% if 1: cumulate the series
% if 2: take the level
% if 3: take the log
% if 4: take the difference
% if 5: take the growth rate
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
i_system = 1;
for i_var = 1:nvar
      if In_or_out(i_var) == 1
            data(:,i_system) = full_data(:,i_var);
            var_names(i_system) = variables(i_var);
            i_system = i_system + 1;
      elseif In_or_out(i_var) == 0
            i_system = i_system;
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