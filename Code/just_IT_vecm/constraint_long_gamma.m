function [c, ceq] = constraint_long_gamma(B,Xi,gam,r,nvar)
% B is the impact matrix
% Xi is the long run matrix from Lutkepohl (2005)
% gam is the impact a la Barsky-Sims
% r is the number of cointegrations
% nvar is the number of variables
% Objective. Set to zero r column of the matrx Xi*B following Lutkepohl (2005)
% The number of elements of Xi*B sets equal to zero is nvar*r
c = [];
ceq = zeros(nvar*r + 1,1); % r*nvar is the number of constraints
obj_matrix = Xi*B*gam;

longrun_constraint = 0;
if longrun_constraint == 1
      ceq = [gam'*gam; obj_matrix] - [1; zeros(nvar*r,1)]; % nvar*r is the number of constraints
else
      ceq = gam'*gam - 1;
end

end