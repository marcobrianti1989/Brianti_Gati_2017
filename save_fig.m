function save_fig(graphics_path, fig_handle, name)

currentFolder = pwd;

h = fig_handle;
set(h,'PaperOrientation','landscape');
set(h,'PaperPosition', [0 0 11 8]); %[1 1 28 19]
current_time = datestr(now);
filename = [name current_time];
cd(graphics_path)
print(filename,'-dpdf')
cd(currentFolder)

