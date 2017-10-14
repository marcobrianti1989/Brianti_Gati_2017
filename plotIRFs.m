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
        title([name, ' on ' , varname])
        hold off
        grid on      
    end
    
    % if you want to print just remove the '%'
    %print(['figure(' num2str(i_shock) ')_RD'] ,'-dpdf','-fillpage')
    
end