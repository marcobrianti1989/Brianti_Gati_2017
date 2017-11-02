   function [c,ceq] = constraint_finite_horizon(gam,R,q) 
       % gam is for this exercise gam3, 6x1
       % R = the effect at some long horizon before multiplying with gam3. (6x6)
       % q = position of relative prices
         c = [];
         ceq = zeros(2,1); % 2 is the number of constraints
         R_rel_prices = R(q,:);
         ceq_row1 = [gam', zeros(1,size(gam,1))];
         ceq_row2 = [zeros(1,size(gam,1)), R_rel_prices];
         ceq = vertcat(ceq_row1, ceq_row2)*vertcat(gam, gam) - [1,0]';  
         
   end