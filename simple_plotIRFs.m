function simple_plotIRFs(IRFs,ub,lb,h,which_shock, print_figs)
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
    
    for i_shock=1:nshocks
        
        % Draw pretty pictures
        figure(i_shock)
        set(gcf,'color','w'); % sets white background color
        set(gcf, 'Position', get(0, 'Screensize')); % sets the figure fullscreen (hopefully)
        for i_var=1:nvar
            subplot(nvar,1,i_var)
            hold on
            x2 = [periods, fliplr(periods)];
            inBetween = [lb(i_var,1:h,which_shock(i_shock)), fliplr(ub(i_var,1:h,which_shock(i_shock)))];
            fill(x2, inBetween, [0.75 0.75 0.75],'LineStyle','none');
            %plot(periods,ub(i_var,1:h,which_shock(i_shock)),'--','linewidth',1.5,'Color','k')
            %plot(periods,lb(i_var,1:h,which_shock(i_shock)),'--','linewidth',1.5,'Color','k')
            plot(zeros(1,h), 'Color','b')
            plot(periods,IRFs(i_var,1:h,which_shock(i_shock)),'linewidth',1,'Color','r')            
            %title([name, ' on ' , varname],'fontsize',18)
            xt = get(gca, 'XTick');
            set(gca, 'FontSize', 14)
            
            hold off
            grid on
        end
        
        % Save figures if you want to
        if strcmp(print_figs, 'yes')
            invoke_export_fig(name, which_ID_strat)
            close all
        end
    end