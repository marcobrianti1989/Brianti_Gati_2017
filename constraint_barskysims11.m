   function [c,ceq] = constraint_barskysims11(gam)
       
         c = [];
         ceq = gam'*gam-1;
         
   end