   function [c,ceq] = constraint_ryan(gam,R,q) 
       % gam is for this exercise gam3, 6x1
       % R = PA, the LR effect before multiplying with gam3. (nvar x nshocks)
       % q = position of relative prices
         c = [];
         ceq = zeros(2,1); % 2 is the number of constraints
         R_rel_prices = squeeze(R(q,:,:));
         if size(R_rel_prices,2) == 1
               R_rel_prices = R_rel_prices';
         end
         ceq = [gam'*gam, (R_rel_prices*gam)'] - [1, zeros(1,size(R_rel_prices,1))];  
         
   end