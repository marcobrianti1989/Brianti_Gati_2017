   function [c,ceq] = constraint_ryan(gam,R,q) 
       % gam is for this exercise gam3, 6x1
       % R = PA, the LR effect before multiplying with gam3. (6x6)
       % q = position of relative prices
         c = [];
         ceq = zeros(2,1); % 2 is the number of constraints
         R_tfp_on_prices = R(q,:);
         ceq_row1 = [gam', zeros(1,6)];
         ceq_row2 = [zeros(1,6), R_tfp_on_prices];
         ceq = vertcat(ceq_row1, ceq_row2)*vertcat(gam, gam) - [1,0]';  
         
   end