function mod=anal_deriv(mod);
% This program copmutes analytical first and second (if approx=2) derivatives of the function f(yp,y,xp,x) with respect to x, y, xp, and yp.  For documentation, see the paper ``Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy Function,'' by Stephanie Schmitt-Grohe and Martin Uribe (JEDC, vol. 28, January 2004, pp. 755-775).
%
%Inputs: f, x, y, xp, yp, approx
%
%Output: Analytical first and second derivatives of f.
%
%If approx is set at a value different from 2, the program delivers the first derivatives of f and sets second derivatives at zero. If approx equals 2, the program returns first and second derivatives of f. The default value of approx is 2.
%Note: This program requires MATLAB's Symbolic Math Toolbox
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001
%
%
%6/24/09: Modified by Ryan Chahrour to work with model object

%Get arguments from the model object
x = mod.X;
y = mod.Y;
xp = mod.XP;
yp = mod.YP;
f = mod.f;
approx = mod.set.approx_deg;


nx = size(x,2);
ny = size(y,2);
nxp = size(xp,2);
nyp = size(yp,2);

n = size(f,2);


%Compute the first and second derivatives of f
fx = jacobian(f,x);
fxp = jacobian(f,xp);
fy = jacobian(f,y);
fyp = jacobian(f,yp);

mod.fx = fx;
mod.fxp = fxp;
mod.fy = fy;
mod.fyp = fyp;

if approx==2

    mod.fypyp = reshape(jacobian(fyp(:),yp),n,nyp,nyp);

    mod.fypy = reshape(jacobian(fyp(:),y),n,nyp,ny);

    mod.fypxp = reshape(jacobian(fyp(:),xp),n,nyp,nxp);

    mod.fypx = reshape(jacobian(fyp(:),x),n,nyp,nx);

    mod.fyyp = reshape(jacobian(fy(:),yp),n,ny,nyp);

    mod.fyy = reshape(jacobian(fy(:),y),n,ny,ny);

    mod.fyxp = reshape(jacobian(fy(:),xp),n,ny,nxp);

    mod.fyx = reshape(jacobian(fy(:),x),n,ny,nx);

    mod.fxpyp = reshape(jacobian(fxp(:),yp),n,nxp,nyp);

    mod.fxpy = reshape(jacobian(fxp(:),y),n,nxp,ny);

    mod.fxpxp = reshape(jacobian(fxp(:),xp),n,nxp,nxp);

    mod.fxpx = reshape(jacobian(fxp(:),x),n,nxp,nx);

    mod.fxyp = reshape(jacobian(fx(:),yp),n,nx,nyp);

    mod.fxy = reshape(jacobian(fx(:),y),n,nx,ny);

    mod.fxxp = reshape(jacobian(fx(:),xp),n,nx,nxp);

    mod.fxx = reshape(jacobian(fx(:),x),n,nx,nx);

% else
% 
%     fypyp=0; fypy=0; fypxp=0; fypx=0; fyyp=0; fyy=0; fyxp=0; fyx=0; fxpyp=0; fxpy=0; fxpxp=0; fxpx=0; fxyp=0; fxy=0; fxxp=0; fxx=0;

end
