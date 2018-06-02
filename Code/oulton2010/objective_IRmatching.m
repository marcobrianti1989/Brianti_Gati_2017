function wlf = objective_IRmatching(param,set,Sy,Sx,H,psi_hat,W)

%Step 1: compute the model solution
[~, fx, fy, fxp, fyp, G, ~, ~]=model_prog(param,set); 
[gx,hx]=gx_hx_alt(fy,fx,fyp,fxp);

nx = size(hx,1);
ny = size(gx,1);

%Step 2: compute the model impulse responses, put y-repsonse on levels
g = Sy*gx;
h = Sx*hx;

%%%%%
% Need to edit the part below:
%%%%%
x0 = 1;  % we only have one shock   % --->>>>> this is specific to BriantiGati2017
[~,ir_IT] = ir(g,h,G*x0,H);  % IT prod level shock

ir_IT(:,:) = cumsum(ir_IT(:,:));


%Step 3: compute the deviations of 
psi_theta = [ir_IT(:)]; % --->>>>> this is specific to BriantiGati2017

%%%%% until here.

dev = psi_hat - psi_theta;


%Step 4: compute the weighted loss function
wlf = 1000*dev'*W*dev;
%disp '----- IR-matching fmincon ------'
