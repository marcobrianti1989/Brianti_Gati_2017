function plot_IRF_lp_unconditional(varlist,IRF_low,IRF_low2,IRF_up,IRF_up2,IRF,H,plot2,n_row,unique)

if plot2 == 1
      %Impulse Response Functions using Local Projection - Figure
      if unique == 1
            nvar     = length(varlist);
            n_col    = ceil(nvar/n_row); %plus one for Vix
            figure('Position',[1 41 1920 963])
            set(gcf,'color','w');
      end      
      for j = 1:length(varlist)
            if unique == 1
                  s = subplot(n_row,n_col,j);
            else
                  figure(j)
                  nvar     = length(varlist);
                  n_col    = ceil(nvar/n_row); %plus one for Vix
                  figure('Position',[1 41 1920 963])
                  set(gcf,'color','w');
            end
            hold on
            plot([0:H-1]',IRF_low(j,:), '--k','linewidth', 1);
            plot([0:H-1]',IRF_up(j,:), '--k','linewidth', 1);
            plot([0:H-1]',IRF_low2(j,:), '--k','linewidth', 2);
            plot([0:H-1]',IRF_up2(j,:), '--k','linewidth', 2);
            plot([0:H-1]',IRF(j,:), '-k', 'linewidth', 3);
            plot([0:H-1]',0*[1:H]',':k');
            set(gca,'TickLabelInterpreter','latex')
            title(varlist{j},'interpreter', 'latex', 'fontsize', 24);
            if unique == 1 && j == 1
                  xlabel('Quarter','interpreter','latex','fontsize',18);
                  ylabel('\% deviation from s.s.','interpreter','latex','fontsize',18);
            elseif unique == 0
                  xlabel('Quarter','interpreter','latex','fontsize',18);
                  ylabel('\% deviation from s.s.','interpreter','latex','fontsize',18);
            end
            if unique == 1
                  set(s,'xlim',[1,H],'ylim', ylim);
            end
            
      end
end