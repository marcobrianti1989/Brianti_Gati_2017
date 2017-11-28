function loss = obj_get_struct_shocks_alt(D_vec, gam,which_shocks)
% Inputs
% D_vec = the vector in which we collect all the elements of D that are not
% gamma, i.e. these elements are the maximand. (nvar*[nvar-size(gam,2)], 1)

nvar = size(gam,1);
D = zeros(nvar, nvar);

for i_var = 1:nvar
      if i_var == which_shocks(1) 
            D(:,i_var) = gam(:,1);
      elseif i_var == which_shocks(2)
            D(:,i_var) = gam(:,2);
      else
            D(:,i_var) = D_vec((i_var*nvar)+1-nvar-(nvar*size(gam,2)):i_var*nvar-(nvar*size(gam,2)));
      end
end

% D(:,1) = D_vec(1:nvar);
% D(:,2) = D_vec(nvar+1:2*nvar);
% D(:,3) = gam(:,1);
% D(:,4) = gam(:,2);
% D(:,5) = D_vec(2*nvar+1:3*nvar);
% D(:,6) = D_vec(3*nvar+1:end);

loss = sum(sum((D'*D - eye(nvar)).^2));
