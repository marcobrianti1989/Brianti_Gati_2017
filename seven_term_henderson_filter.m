function [dt, seasonality] = seven_term_henderson_filter(y)

wgh      = [-0.059 0.059 0.294 0.412 0.294 0.059 -0.059];
T        = length(y);

c        = 4;
dt(1)    = (y(1)*wgh(c) + y(2)*wgh(c+1) + y(3)*wgh(c+2))...
      /(wgh(c) + wgh(c+1) + wgh(c+2));
dt(2)    = (y(1)*wgh(c-1) + y(2)*wgh(c) + y(3)*wgh(c+1) + y(4)*wgh(c+2))...
      /(wgh(c-1) + wgh(c) + wgh(c+1) + wgh(c+2));
dt(3)    = (y(1)*wgh(c-2) + y(2)*wgh(c-1) + y(3)*wgh(c) + y(4)*wgh(c+1) + y(5)*wgh(c+2))...
      /(wgh(c-1) + wgh(c) + wgh(c+1) + wgh(c+2));
for i = 4:T-3
   dt(i) = wgh(c-3)*y(i-3) + wgh(c-2)*y(i-2) + wgh(c-1)*y(i-1) ...
         + wgh(c)*y(i) + wgh(c+1)*y(i+1) + wgh(c+2)*y(i+2) + wgh(c+3)*y(i+3);
end
dt(T-2)    = (y(T-4)*wgh(c-2) + y(T-3)*wgh(c-1) + y(T-2)*wgh(c) + y(T-1)*wgh(c+1) + y(T)*wgh(c+2))...
      /(wgh(c-2) + wgh(c-2) + wgh(c) + wgh(c+1) + wgh(c+2));
dt(T-1)    = (y(T-3)*wgh(c-2) + y(T-2)*wgh(c-1) + y(T-1)*wgh(c) + y(T)*wgh(c+1))...
      /(wgh(c-2) + wgh(c-2) + wgh(c) + wgh(c+1));
dt(T)      = (y(T-2)*wgh(c-2) + y(T-1)*wgh(c-1) + y(T)*wgh(c))...
      /(wgh(c-2) + wgh(c-2) + wgh(c));

seasonality = y - dt;

end