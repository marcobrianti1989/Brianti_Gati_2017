function sim_data = simulate_model(gx,hx,x0,T)
% builds on the structure of ir.m of Ryan to generate simulated time series
% from the model solution gx,hx.
% x0 should be eta*(shock vector of zeros and 1 for the shock in question).

if nargin < 4
    T = 100;
end

x0=x0(:);
pd=length(x0);
MX=[gx;eye(pd)]; % MX is a vector of [gx(:); hx(:)] which in each iteration gets multiplied by the impulse of that period.
sim_data=[]; % (T, nvar) with the jumps first and then the states.
% x=x0;
for t=1:T
    x = x0*randn(1); % only modification vis-a-vis ir.m : every period there is a shock of rand size.
    sim_data(t,:)=(MX*x)';
    x = hx * x;
end
ny = size(gx,1);
sim_y = sim_data(:,1:ny); %jumps
sim_x = sim_data(:,ny+1:end); %states