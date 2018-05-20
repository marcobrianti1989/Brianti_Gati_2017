%
% mom_tab - make the basic moment tables to screen
%
% usage:
% 
% mom_tab(gx, hx, varshock, var_idx, var_names)
%
% where
%
% gx,hx,varshock are as in mom.m
% var_idx = indexies of the relevant variables in Y (first value should index output)
% var_names = names of the variables

function [sig ar_rho rho] = mom_tab(gx, hx, varshock, var_idx, var_names)




[SigY0, SigX0] = mom(gx, hx, varshock);
[SigY1, SigX1] = mom(gx, hx, varshock, -1);


yidx = var_idx(1);
tit_str = [];
for j = 1:length(var_idx)
    idx = var_idx(j);
    sig(j) = sqrt(SigY0(idx,idx));
    rho(j) = SigY0(idx,yidx)/sqrt(SigY0(idx,idx)*SigY0(yidx,yidx));
    ar_rho(j) = SigY1(idx,idx)/(sig(j)^2);
    tit_str = [tit_str, var_names{j}, '\t'];
end


    
disp('Standard Deviations')
disp(sprintf(tit_str))
disp(sprintf('%1.3f\t', 100*sig))
disp(' ');


disp('Stddev(X)/Stddev(Y)')
disp(sprintf(tit_str))
disp(sprintf('%1.3f\t', sig/(sig(1))))
disp(' ');


disp('Auto-correlations')
disp(sprintf(tit_str))
disp(sprintf('%1.3f\t', ar_rho))
disp(' ');

disp('Correlation w/ Y')
disp(sprintf(tit_str))
disp(sprintf('%1.3f\t', rho))
disp(' ');
