function plot_matchedIRFs(IRFs_VAR, IRFs_theory,h,which_shock,shocknames, varnames, print_figs, base_path, comment)
% Inputs
% IRFs_VAR    = IRFs from VAR
% IRFs_theory = IRFs from IR-matching exercise
% h = horizon
% which_shock = a vector or scalar specifying which are the shocks
% shocknames  = cell vector of names of shocks
% varnames    = cell vector of names of variables
% print_figs  = 'yes' or 'no'
% base_path   = should really be the base path; the code automatically puts
% figs in the Figures folder
% comment     = put in something so that when saved, something
% distinguishes the figures. If empty, '_' will be put in.
% Code written by Marco Brianti & Laura Gáti, 24 June 2018.

nvar = size(IRFs_VAR,1);
nshocks = size(which_shock,2);
periods = 1:h;

%Ylim
min_y_lim = [-0.5    -0.1     -0.01     0        0        -0.012];
max_y_lim = [2.5      0.8      0.03     0.012    0.012     0.006];

% Draw pretty pictures
for i_shock=1:nshocks
    for i_var=1:nvar
        figure(i_shock*nvar - nvar + i_var)
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
        varname = varnames{i_var};
        name = shocknames{which_shock(i_shock)};
        hold on
        plot(periods,IRFs_VAR(i_var,1:h,which_shock(i_shock)),'linewidth',2,'Color','k')
        plot(periods,IRFs_theory(i_var,1:h,which_shock(i_shock)),'--','linewidth',2,'Color','r')
        plot(zeros(1,h), 'Color','b')
        LEG = legend('Empirical', 'Model');
        LEG.FontSize = 24;
        legend boxoff
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 22)
        title([name, ' on ' , varname],'fontsize',72)
        %   ylim([min_y_lim(i_var) max_y_lim(i_var)]);
        hold off
        grid on
        
        % Save figures if you want to
        if strcmp(print_figs, 'yes')
            if nargin < 8
                comment = '_';
            end
            use_current_time = 'no';
            invoke_export_fig([name, ' on ' , varname], comment,use_current_time, base_path)
            close all
            pause(0.5)
        end
    end
end


end