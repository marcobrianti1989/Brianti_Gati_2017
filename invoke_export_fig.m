function invoke_export_fig(figname, some_comment)
% This file invokes the export_fig command which requires the figure to be saved to be open
% Inputs:
% figname      = what you want the figure to be called
% some_comment = something that will show up in the figure title, e.g. ID strategy
% Outputs: a figure saved as a pdf with the name
% fig_figname_some_comment_time.pdf


current_time = datestr(now);
ct1 = strrep(current_time,':','_'); % substitute : with _
ct2 = strrep(ct1,' ','_'); % substitute spaces with _

% Make sure there are no spaces or : in any of the strings
figname = strrep(figname,':','_');
figname = strrep(figname,' ','_');
some_comment = strrep(some_comment,':','_');
some_comment = strrep(some_comment,' ','_');

final_name = ['fig_' figname '_' some_comment '_' ct2 '.pdf'];

base_path = pwd;
if exist([base_path '\Figures'], 'dir')
    cd([base_path '\Figures']) %for Microsoft
else
    cd([base_path '/Figures']) %for Mac
end

if exist([base_path '\Export_Fig'], 'dir')
    addpath([base_path '\Export_Fig']) %for Microsoft
else
    addpath([base_path '/Export_Fig']) %for Mac
end

export_fig(final_name)

cd(base_path)


end