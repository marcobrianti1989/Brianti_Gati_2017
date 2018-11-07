function export_fig_IRF_lp_unconditional(export_fig2)

if export_fig2 == 1
      
      % Create the correct path
      base_path = pwd;
      warning off
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
      warning on
      export_fig(['DifferentInflationControls_logLevels_8_30_2018.pdf'])
      close all
      cd(base_path) %back to the original path

end



end