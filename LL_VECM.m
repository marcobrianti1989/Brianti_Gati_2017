function obj = LL_VECM(T,B,sigma)    
    %obj = T/2*log(det(B)^2) + sum(sum(abs(T/2*trace((B')^(-1)*B^(-1)*sigma))));
    obj = sum(sum((B*B'-sigma).^2));
end