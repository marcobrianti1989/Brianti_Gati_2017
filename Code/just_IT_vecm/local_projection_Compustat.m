clear
close all

%*************************************************************************%
%                                                                         %
%                         Local Projections                               %
%                                                                         %
%*************************************************************************%

%Data Reading and Transformation
filename = 'Quarterly';
sheet    = 'Quarterly Data';
range    = 'B1:CJ275';
warning off
[data, varnames] = read_data2(filename, sheet, range);
warning on
for i = 1:size(data,2)
      eval([varnames{i} ' = data(:,i);'])
end

%Read dataset_PC for PC analysis
filename_PC                                = 'Dataset_test_PC';
sheet_PC                                   = 'Quarterly';
range_PC                                   = 'B2:DC288';
do_truncation_PC                           = 1; %Do truncate data.
[dataset_PC, var_names_PC]                 = read_data2(filename_PC,sheet_PC,range_PC);
dataset_PC                                 = real(dataset_PC);
date_start_PC                              = dataset_PC(1,1);
dataset_PC                                 = dataset_PC(:,2:end); %Removing time before PC analysis
Zscore                                     = 1; %Standardize data before taking PC
PC                                         = get_principal_components(dataset_PC);
mpc                                        = 4;
pckk                                       = PC(:,1:mpc);

% Create Var List
varlist          = {'RD2Assets', 'RD2Investment','RD2Equity',...
      'Soft2Assets','Soft2Investment','Soft2Equity'};
numberRD2A        = strmatch('RD2Assets', varlist);
numberRD2I        = strmatch('RD2Investment', varlist);
numberRD2E        = strmatch('RD2Equity', varlist);
numberS2A         = strmatch('Soft2Assets', varlist);
numberS2I         = strmatch('Soft2Investment', varlist);
numberS2E         = strmatch('Soft2Equity', varlist);

for i = 1:length(varlist)
      dep_var(:,i) = eval(varlist{i});
end

% Loading ICT shock
load ssICT_time
timeICTss = ssICT_time(:,1);
ICTss     = ssICT_time(:,2);

% Align datasets
if timeICTss(1) < Time(1) && timeICTss(end) >= Time(end)
      disp(['Alignment is Case 1'])
      loc_start       = find(timeICTss(:) == Time(1));
      loc_end         = find(timeICTss(:) == Time(end));
      ICTss           = ICTss(loc_start:loc_end);
      pckk            = pckk(loc_start:loc_end,:);
elseif timeICTss(1) > Time(1) && timeICTss(end) <= Time(end)
      disp(['Alignment is Case 2'])
      loc_start       = find(Time(:) == timeICTss(1));
      loc_end         = find(Time(:) == timeICTss(end));
      dep_var         = dep_var(loc_start:loc_end);
      pckk            = pckk(loc_start:loc_end,:);
elseif timeICTss(1) < Time(1) && timeICTss(end) <= Time(end)
      disp(['Alignment is Case 3'])
      loc_start       = find(timeICTss(:) == Time(1));
      loc_end         = find(Time(:) == timeICTss(end));
      ICTss           = ICTss(loc_start:end);
      dep_var         = dep_var(1:loc_end);
      pckk            = pckk(loc_start:loc_end,:);
elseif timeICTss(1) > Time(1) && timeICTss(end) >= Time(end)
      disp(['Alignment is Case 4'])
      loc_start       = find(Time(:) == timeICTss(1));
      loc_end         = find(timeICTss(:) == Time(end));
      ICTss           = ICTss(1:loc_end);
      dep_var         = dep_var(loc_start:end);
      pckk            = pckk(loc_start:loc_end,:);
end

%numberInflation  = strmatch('Inflation', varlist);
lags             = 2;
H                = 20; %irfs horizon

% Standardize Ztilde to get one std dev shock
ss  = ICTss/std(ICTss);

% Matrix of dependen variables - All the variables are in log levels
control_pop = 0; % Divide GDP, Cons, Hours, Investment over population
for i = 1:length(varlist)
      dep_var(:,i) = eval(varlist{i});
      if control_pop == 1
            if i == numberGDP || i == numberC || i == numberHours || i == numberInv || i == numberInvent || i == numberInvent
                  dep_var(:,i) = dep_var(:,i) - Population;
            end
      end
end

% Set up the typology of transformation
logdifferences = 0;
if logdifferences == 1
      dep_var = [nan(1,size(dep_var,2)); diff(dep_var)];
end

for kk = 1:size(dep_var,2)
      % Define inputs for local_projection
      depvarkk                    = dep_var(:,kk);
      % Run local_projection
      [IR{kk},res{kk},Rsquared{kk},BL{kk},tuple{kk},VarY{kk}] = ...
            local_projection(depvarkk,pckk,ss,lags,H);
      if logdifferences == 0
            IRF(kk,:) = IR{kk};
      else
            IRF(kk,:) = cumsum(IR{kk});
      end
      % Build a table for the Variance Explained by Ztilde - Following  Stock,
      % Watson (2018) - The Economic Journal, page 928 Eq. (15)
      VarY_ih = VarY{kk};
      for ih = 1:H
            VarYY    = VarY_ih(ih);
            VarExplained(kk,ih) = sum(IRF(kk,1:ih).^2)/VarYY;
      end
      % Initiate bootstrap
      nsimul         = 2000;
      tuplekk        = tuple{kk};
      for hh = 1:H
            tuplekkhh = tuplekk{hh}; % Fix a specific horizon
            Y                             = tuplekkhh(:,1);
            X                             = tuplekkhh(:,2:end);
            XControl                      = tuplekkhh(:,3:end);
            [Yboot, Xboot]                = bb_bootstrap_LP(Y,X,nsimul,lags);
            [YbootC, XbootC]              = bb_bootstrap_LP(Y,XControl,nsimul,lags);
            for isimul = 1:nsimul
                  B                       = Xboot(:,:,isimul)'*Xboot(:,:,isimul)\...
                        (Xboot(:,:,isimul)'*Yboot(:,isimul));
                  BC                      = XbootC(:,:,isimul)'*XbootC(:,:,isimul)\...
                        (XbootC(:,:,isimul)'*YbootC(:,isimul));
                  IRF_boot(kk,hh,isimul)  = B(1);
                  VarYBoot(kk,hh,isimul)  = var(YbootC(:,isimul) - XbootC(:,:,isimul)*BC);
            end
      end
end

% Select upper and lower bands
for kk = 1:size(dep_var,2)
      IRF_bootkk = IRF_boot(kk,:,:);
      VarYbootkk = VarYBoot(kk,:,:);
      if logdifferences == 0
            IRF_boot(kk,:,:)  = IRF_bootkk;
            VarY_boot(kk,:,:) = VarYbootkk;
      else
            IRF_boot(kk,:,:)  = cumsum(IRF_bootkk,2);
            VarY_boot(kk,:,:) = cumsum(VarYbootkk,2);
      end
end
IRF_boot         = sort(IRF_boot,3);
VarY_boot        = sort(VarY_boot,3);
sig              = 0.05;
sig2             = 0.16;
up_bound         = floor(nsimul*sig); % the upper percentile of bootstrapped responses for CI
up_bound2        = floor(nsimul*sig2); % the upper percentile of bootstrapped responses for CI
low_bound        = ceil(nsimul*(1-sig)); % the lower percentile of bootstrapped responses for CI
low_bound2       = ceil(nsimul*(1-sig2)); % the lower percentile of bootstrapped responses for CI
IRF_up           = IRF_boot(:,:,up_bound);
VarY_up          = VarY_boot(:,:,up_bound);
IRF_up2          = IRF_boot(:,:,up_bound2);
VarY_up2         = VarY_boot(:,:,up_bound2);
IRF_low          = IRF_boot(:,:,low_bound);
VarY_low         = VarY_boot(:,:,low_bound);
IRF_low2         = IRF_boot(:,:,low_bound2);
VarY_low2        = VarY_boot(:,:,low_bound2);

% Confidence Intervals for Variance Explained
for kk = 1:size(dep_var,2)
      VarYup   = VarY_up(kk,:);
      VarYup2  = VarY_up2(kk,:);
      VarYlow  = VarY_low(kk,:);
      VarYlow2 = VarY_low2(kk,:);
      for ih = 1:H
            VarYYup   = VarYup(ih);
            VarYYup2  = VarYup2(ih);
            VarYYlow  = VarYlow(ih);
            VarYYlow2 = VarYlow2(ih);
            VarExplainedup(kk,ih)   = sum(IRF_up(kk,1:ih).^2)/VarYYup;
            VarExplainedlow(kk,ih)  = sum(IRF_low(kk,1:ih).^2)/VarYYlow;
            VarExplainedup2(kk,ih)  = sum(IRF_up2(kk,1:ih).^2)/VarYYup2;
            VarExplainedlow2(kk,ih) = sum(IRF_low2(kk,1:ih).^2)/VarYYlow2;
      end
end

%Show the graph of IRF - Figure(2)
plot2    = 1; % if plot2 = 1, figure will be displayed
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,plot2,n_row,unique)

% Print figure authomatically if "export_figure1 = 1"
if plot2 == 1
      export_fig2 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_IRF_lp_unconditional(export_fig2)
end
asd
%Show the variance Explained - Figure(3)
plot3    = 1; % if plot2 = 1, figure will be displayed
n_row    = 3; % how many row in the figure
unique   = 1; % if unique = 1 plot IRFs together, if = 1 plot each IRF separately
plot_IRF_lp_unconditional(varlist,VarExplained,VarExplained,VarExplained,...
      VarExplained,VarExplained,H,plot3,n_row,unique)

% Print figure authomatically if "export_figure1 = 1"
if plot3 == 1
      export_fig3 = 0; % if export_fig1 = 1, figure will be saved
      export_fig_IRF_lp_unconditional(export_fig3)
end
