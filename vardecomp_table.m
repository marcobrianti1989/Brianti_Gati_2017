function [vardec_table] = vardecomp_table(vardec,which_shock,varnames,names)

nvar = size(vardec,1);
nshocks = size(which_shock,2);
vardec_table = cell(nvar+1,nshocks+1);
vardec_table(1,2:end) = names;    
vardec_table(2:end,1) = varnames';
vardec_table(2:end,2:end) = num2cell(vardec(:,which_shock))

end