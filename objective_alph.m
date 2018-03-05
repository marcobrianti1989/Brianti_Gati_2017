function obj_alph = objective_alph(Pi,r,alph_bet) 
% This function is simpli Pi - alp*bet' from the cointegration matrix, the
% cointegrating verctors bet(i) and the loading matrix alph.

alph = alph_bet(:,r);
bet  = alph_bet(:,r+1:end);

%Objective
obj_alph1  = Pi - alph*bet';
obj_alph   = sum(sum(obj_alph1.^2));

end