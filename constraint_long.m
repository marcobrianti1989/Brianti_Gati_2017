function [c, ceq] = constraint_long(B,Xi,r,nvar)
% B is the impact matrix
% Xi is the long run matrix from Lutkepohl (2005)
% r is the number of cointegrations
% nvar is the number of variables
% Objective. Set to zero r column of the matrx Xi*B following Lutkepohl (2005)
% The number of elements of Xi*B sets equal to zero is nvar*r
c = [];
obj_matrix = Xi*B;
obj_matrix(:,nvar-r:end)
ceq = vertcat(obj_matrix(:,nvar-r:end)); % nvar*r is the number of constraints

end