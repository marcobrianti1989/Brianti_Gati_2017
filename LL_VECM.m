function obj = LL_VECM(T,B,sigma)    
    obj = sum(sum(abs(T/2*log(B.^2)+T/2*trace((B')^(-1)*B^(-1)*sigma))));   
end