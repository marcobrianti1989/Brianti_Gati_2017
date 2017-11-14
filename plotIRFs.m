function plotIRFs(IRFs,ub,lb,h,which_shock, names, varnames)
    % h = IR horizon for figures
    % ub and lb are the bootstrapped CI
    % names = cell vector of shock names
    % varnames = cell vector of variable names
    
    nvar = size(IRFs,1);
    nshocks = size(which_shock,2);
    periods = 1:h;
    
    for i_shock=1:nshocks
        name = names{i_shock};
        
        % Draw pretty pictures
        figure(i_shock)
        for i_var=1:nvar
            varname = varnames{i_var};
            subplot(nvar,1,i_var)
            hold on
            plot(periods,IRFs(i_var,1:h,which_shock(i_shock)),'linewidth',1.5,'Color','r')
            plot(periods,ub(i_var,1:h,which_shock(i_shock)),'--','linewidth',1,'Color','k')
            plot(periods,lb(i_var,1:h,which_shock(i_shock)),'--','linewidth',1,'Color','k')
            plot(zeros(1,h), 'Color','b')
            title([name, ' on ' , varname])
            hold off
            grid on
        end
        
        % if you want to print just remove the '%'
        figname = ['figure_' , datestr(now, 'yyyy-mm-dd HH:MM:SS')];
        %     print(figname ,'-dpdf') %,'-fillpage')
        %     pause(2)
        
        % print(['figure_Nov14_', names{i_shock}] ,'-dpdf') %,'-fillpage')
        
        
    end