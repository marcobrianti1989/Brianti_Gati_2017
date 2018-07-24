function known_ss_values = get_ss_ratios_from_data

current_dir = pwd;
cd ../.. % go up 2 levels
base_path = pwd;
cd(current_dir)
addpath(base_path)

%Data Reading and Transformation (I'm taking growth rates of nonstationary
%variables)
filename = 'SS_averages';
sheet    = 'data';
range    = 'C1:L118';
[data, varnames] = read_data3(filename, sheet, range, base_path);
varnames

% Take averages
rc = mean(data(:,1));
p  = mean(data(:,2));
c  = mean(data(:,3));
ic = mean(data(:,4));
it = mean(data(:,5));
w  = mean(data(:,6));
ri = mean(data(:,7));
h  = mean(data(:,8));
kc = mean(data(:,9));

known_ss_values = [ic it c w rc ri p h kc];