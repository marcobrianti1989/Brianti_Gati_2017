  function [c,ceq] = constraint_orthogonality_ryan(gam)
       
         c = [];
         ceq = gam'*gam-eye(size(gam,2));
         
   end