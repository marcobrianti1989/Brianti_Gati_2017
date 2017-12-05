function plot_single_IRFs(IRFs,ub,lb,h,which_shock,names, varnames, which_ID_strat, print_figs)
% h = IR horizon for figures
% ub and lb are the bootstrapped CI
% names = cell vector of shock names
% varnames = cell vector of variable names
% which_ID_strat = a string describing which identification strategy we used
% print_figs = 'yes' --> saves the figures; else if 'no' -- > just
% shows the figures.

nvar = size(IRFs,1);
nshocks = size(which_shock,2);
periods = 1:h;

% Draw pretty pictures
for i_shock=1:nshocks
      for i_var=1:nvar
            figure(i_shock*nvar - nvar + i_var)
            set(gcf,'color','w'); % sets white background color
            set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
            varname = varnames{i_var};
            name = names{i_shock};
            hold on
            x2 = [periods, fliplr(periods)];
            inBetween = [lb(i_var,1:h,which_shock(i_shock)), ...
                  fliplr(ub(i_var,1:h,which_shock(i_shock)))];
            fill(x2, inBetween, [0.75 0.75 0.75],'LineStyle','none');
            plot(zeros(1,h), 'Color','b')
            plot(periods,IRFs(i_var,1:h,which_shock(i_shock)),'linewidth',1,'Color','r')
            %title([varname],'fontsize',18)
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', 14)
            title([name, ' on ' , varname],'fontsize',60)
            hold off
            grid on
            
            % Save figures if you want to
            if strcmp(print_figs, 'yes')
                  invoke_export_fig([name, ' on ' , varname], which_ID_strat)
                  close all
                  pause(0.5)
            end
      end
end


end