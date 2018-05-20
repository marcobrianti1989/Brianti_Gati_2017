function wlf = objective_IRmatching(param,set,S,H,psi_hat,W)

%Step 1: compute the model solution
[~, fx, fy, fxp, fyp, G, ~, ~]=model_prog_IRmatching_spillover_news(param,set); % --->>>>> this is specific to BriantiGati2017
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);


%Step 2: compute the model impulse responses, put y-repsonse on levels
g = S*gx;

%%%%%
% Need to edit the part below:
%%%%%
[~,ir_tfp] = ir(g,hx,G*[1,0]',H);  %TFP Shock
[~,ir_mon] = ir(g,hx,G*[0,1]',H);  %Monetary Shochk

ir_tfp(:,1) = cumsum(ir_tfp(:,1));
ir_mon(:,1) = cumsum(ir_mon(:,1));


%Step 3: compute the deviations of 
psi_theta = [ir_tfp(:); ir_mon(:)];

%%%%% until here.

dev = psi_hat - psi_theta;


%Step 4: compute the weighted loss function
wlf = 1000*dev'*W*dev;
