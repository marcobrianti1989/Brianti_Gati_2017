function [structural_shocks, IR] = get_structural_shocks_general(A,gamma,resid,which_shocks)

D = [gamma null(gamma')]; % building D, the only orthogonal matrix
%D = [D_step1(:,1:which_shocks(1)-1) gamma D_step1(:,which_shocks(end)+1:end)]; % rearranging so the identified shocks are still in the right place
IR = A*D; 
s = IR^(-1)*resid';
s = s'; % to make it (T,nshocks_recovered)
structural_shocks = s(:,1:length(which_shocks));
end