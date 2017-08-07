%Model 1 steady state (computed analytically)

function [ss] = model_1_ss()
    
    param
    
    %Use closed form expressions for the ss values.
    r = 1/bet - 1 + delt;
    y_k = r/((1-gamm)*alph);
    y_m = nu/gamm;
    wl_y = (1-gamm)*(1-alph);
    v_y = 1/(1-phi*bet)*(nu-1)/y_m;
    lamby = 1/(phi*bet*v_y);
    s_y = (lambx*phi*bet*v_y*y_k^rho)^(1/(1-rho));
    z = lamby/(1-phi)*s_y;
    c_y = - s_y + wl_y + (1-bet)/bet/y_k + (nu - 1)/y_m;
    l = (chi/c_y*wl_y)^(1/thet);
    y = (y_k^(-alph*(1-gamm))*l^((1-alph)*(1-gamm))*z^(nu*gamm)*y_m^(-gamm))^(1/((1-alph)*(1-gamm)));
    c = c_y*y;
    s = s_y*y;
    lamb = lamby/y;
    v = v_y*y;
    w = wl_y*y/l;
    m = y*gamm/nu;
    k = y/y_k;
    pi = (nu-1)*m;
    a = 1;
    
    
    
    %Put the ss values in a vector consistent with Y and X vectors in model.m
    yy  = [y c l s w r m lamb pi v];
    xx  = [k z];
    ss  = [yy xx];
    
end