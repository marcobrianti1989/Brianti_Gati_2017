function [c,ceq] = constraint_FEVmax_sr(gam,A,q)
    % gam is for this exercise gam3, 6x1
    % A = impact matrix coming from Choleski (nvar,nvar)
    % q = position of Mich Index
    
    nvar = size(A,1);
    
    c = [];
    ceq = zeros(4,2); % 3 is the number of constraints: 2x2 for orthognality, 1 for impact of IT on Mich 
    %and a fake fourth one to agree with matrix dimension
    gam = [gam(1:nvar,1) gam(nvar+1:end,1)]; %Unvectorize gam
    gam4 = [gam(:,2) zeros(nvar,1)]; %pad with zeros for matrix dimension alignment
    Mich_indx = [A(q,:); zeros(1,nvar)]; %Impacts on Mich from cholesky matrix
    ceq_row12 = [gam', zeros(2,nvar)]; 
    ceq_row34 = [zeros(2,nvar), Mich_indx];
    ceq = vertcat(ceq_row12, ceq_row34)*vertcat(gam, gam4) - [eye(2), zeros(2,2)]';
    
end