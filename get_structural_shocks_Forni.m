function [s, IR] = get_structural_shocks_Forni(A,gamma,resid)

D_step1 = [gamma null(gamma')]; % building D, the only orthogonal matrix
D = [D_step1(:,3) gamma D_step1(:,4:6)]; % rearranging so the identified shocks are still in the right place
IR = A*D; 
s = IR^(-1)*resid';
s = s'; % to make it (T,nshocks_recovered)

end