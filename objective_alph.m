function obj_alph = objective_alph(Pi,bet,alph_vec) 
% This function is simpli Pi - alp*bet' from the cointegration matrix, the
% cointegrating verctors bet(i) and the loading matrix alph.

%Objective
obj_alph1  = Pi - alph_vec*bet';
obj_alph   = sum(sum(obj_alph1.^2));

end