% MODEL_FUNC - Convert a mod object into a matlab file that evaluates the
% first order derivatives of a model without reference to the symbolic
% toolbox.
%
% Usage:
%
% model_func(mod_obj,trans,lb,ub)
%
% where
%
% mod_obj = a data structure which contains the following fields...
%
%
% trans (optional) = an vector of lenth N(parameters) that specifies how to transform the paramters: 0-No Transformation,
%                    1-Transform to restrict the paramters to be > 0, 2-Transform to restrict between lb and ub
%
% lb (optional)   = lower bounds for variables which have transformation 2
% ub (optional)   = upper bounds for variables which have transformation 2
%
% Code by Ryan Chahrour, Boston College, 2012

function model_func(model, varargin)

%**************************************************************************
%This section of the code generates matlab programs for transforming/untransforming the
%parameter vector so that parameters can easily be bound, either above zero
%or between a and b
%**************************************************************************
if ~isempty(varargin)
    f1 = fopen('param_trans.m', 'w');
    str = ['function param = param_trans(param)'];
    fprintf(f1, '%s\n', str);fprintf(f1, '%s\n', '');
    
    
    f2 = fopen('param_untrans.m', 'w');
    str = ['function param = param_untrans(param)'];
    fprintf(f2, '%s\n', str);fprintf(f2, '%s\n', '');
    trans = varargin{1};
    
    %Transform the parameters if requested
    for j =1:length(trans)
        pstr = ['param(',num2str(j) ')'];
        switch trans(j)
            case  1
                %Bound above 0
                str = [pstr '=exp(' pstr ');'];
                fprintf(f1, '%s\n', str);
                
                str = [pstr '=log(' pstr ');'];
                fprintf(f2, '%s\n', str);
            case 2
                %Bound between a and b
                a = num2str(varargin{2}(j));
                b = num2str(varargin{3}(j));
                str = [pstr '=' a '+(' b '-' a ')*exp(' pstr ')/(1+exp(' pstr,'));'];
                fprintf(f1, '%s\n', str);
                
                str = [pstr '=real(log((' pstr '-' a ')/(' b '-' pstr,')));'];
                fprintf(f2, '%s\n', str);
        end
    end
    fclose(f1);
    fclose(f2);
end

f = fopen(model.fname, 'w');

fp = fopen('param_unpacker.m', 'w'); %Script to unpack parameters

%The function call will take the parameters and settings as arguments
str = 'function [f fx fy fxp fyp eta R set dgam_dtheta deta_dtheta dR_dtheta xlag ylag] = model_prog(param, set)';
fprintf(f, '%s\n', str);
fprintf(f, '%s\n', '');

%Put in the transformation
if ~isempty(varargin)
    fprintf(f,'%s\n','param = param_trans(param);');
end

%Assign parameter values to named variables
str_sv{1} = ['function ' str(1:end-1)];
str = '%Assign parameter values to named variables.';
str_sv{2}= str;
fprintf(f, '%s\n', str); fprintf(fp, '%s\n', str);
for j = 1:length(model.PARAM)
    str = [char(model.PARAM(j)) ' = param(' num2str(j), ');'];
    str_sv{j+2} = str;
    fprintf(f, '%s\n',str); fprintf(fp, '%s\n', str);
end
fprintf(f, '%s\n', ''); fprintf(fp, '%s\n','');


%Assign set values to named variables
str = '%Assign set values to named variables.';
str_sv{j+3} = str;
fprintf(f, '%s\n', str); fprintf(fp, '%s\n', str);
for j = 1:length(model.SET)-2
    str = [char(model.SET(j)) ' = set(' num2str(j), ');'];
    str_sv{j+length(model.PARAM)+4} = str;
    fprintf(f, '%s\n',str); fprintf(fp, '%s\n', str);
end
fprintf(f, '%s\n', ''); fprintf(fp, '%s\n', '');


if isfield(model, 'xtra')
    
    %Assign values to matrix parameters
    str = '%Assign vectors values to named vector parameters.';
%    str_sv{j+3} = str;
    fprintf(f, '%s\n', str); fprintf(fp, '%s\n', str);
    
    xlist = fieldnames(model.xtra);
    for j = 1:length(xlist)
        stmp = eval(['symmat_print(model.xsym.' xlist{j} ');']);
        str = [xlist{j} ' = ' stmp  ';'];
        fprintf(f, '%s\n',str); fprintf(fp, '%s\n', str);
    end
    
end

fclose(fp);
%***************************
%Compute SS variables: using external code if needed
%***************************
if any(strcmp('ss_call',fieldnames(model)))
    sscode = fopen(model.ss_call, 'r');
    l = fgetl(sscode);
    j = 1;
    while ~strcmp('%BEGIN_EXTRACT_HERE', l) && j < 300
        l = fgetl(sscode);
        j = j+1;
    end
    
    fprintf(f, '%s\n', l);
    while ~strcmp('%END_EXTRACT_HERE', l) && j < 600
        l = fgetl(sscode);
        j = j+1;
        fprintf(f, '%s\n', l);
    end
    
    fclose(sscode);
    
    
    str = '%Compute Steady State';
    fprintf(f, '%s\n', str);
    for j = 1:length(model.X)
        str = [char(model.X(j)) '= Xss(' num2str(j) ');'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.Y)
        str = [char(model.Y(j)) '= Yss(' num2str(j) ');'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.XP)
        str = [char(model.XP(j)) '= Xss(' num2str(j) ');'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.YP)
        str = [char(model.YP(j)) '= Yss(' num2str(j) ');'];
        fprintf(f, '%s\n', str);
    end
    fprintf(f, '%s\n', '');
    
    
    
else
    str = '%Compute Steady State, testing';
    fprintf(f, '%s\n', str);
    for j = 1:length(model.X)
        str = [char(model.X(j)) '=' char(model.Xss(j)) ';'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.Y)
        str = [char(model.Y(j)) '=' char(model.Yss(j)) ';'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.XP)
        str = [char(model.XP(j)) '=' char(model.X(j)) ';'];
        fprintf(f, '%s\n', str);
    end
    
    for j = 1:length(model.YP)
        str = [char(model.YP(j)) '=' char(model.Y(j)) ';'];
        fprintf(f, '%s\n', str);
    end
    fprintf(f, '%s\n', '');
end

%***************************
%F
%***************************
str = '%Evaluate F.';
fprintf(f, '%s\n', str);

model.f = subs(model.f, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
    [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]);
% L: taking out the log for now and ask Ryan about this!
% model.f = subs(model.f, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
%     log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));

str = ['f = ' symmat_print(model.f) ';'];
fprintf(f, '%s\n', str);

%***************************
%FX
%***************************
str = '%Evaluate derivative expressions.';
fprintf(f, '%s\n', str);
model.fx = subs(model.fx, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
    [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]);
% L: taking out the log for now and ask Ryan about this!
% model.fx = subs(model.fx, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
%     log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));
str = ['fx = ' symmat_print(model.fx) ';'];
fprintf(f, '%s\n', str);


%***************************
%FY
%***************************
model.fy = subs(model.fy, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
    [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]);
% L: taking out the log for now and ask Ryan about this!
% model.fy = subs(model.fy, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
%     log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));
str = ['fy = ' symmat_print(model.fy) ';'];
fprintf(f, '%s\n', str);

%***************************
%FXP
%***************************
model.fxp = subs(model.fxp, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
    [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]);
% L: taking out the log for now and ask Ryan about this!
% model.fxp = subs(model.fxp, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
%     log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));
str = ['fxp = ' symmat_print(model.fxp) ';'];
fprintf(f, '%s\n', str);


%***************************
%FYP
%***************************
model.fyp = subs(model.fyp, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
    [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]);
% L: taking out the log for now and ask Ryan about this!
% model.fyp = subs(model.fyp, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
%     log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));
str = ['fyp = ' symmat_print(model.fyp) ';'];
fprintf(f, '%s\n', str);
fprintf(f, '%s\n', '');


%***************************
%Standard Error Matrices
%***************************
str = ['eta = ' symmat_print(model.shck) ';'];
fprintf(f, '%s\n', str);

str = ['R = ' symmat_print(model.me) ';'];
fprintf(f, '%s\n', str);

%********************************
%INCLUDE ANALYTICAL DERIVATIVES?
%********************************
if model.adiff
    
    %***************************
    %dgam_dtheta (model written at GAM*A = 0)
    %***************************
    model.dgam_dtheta = subs(model.dgam_dtheta, [model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)],...
        log([model.X(model.xlog), model.Y(model.ylog),model.XP(model.xlog), model.YP(model.ylog)]));
    str = ['dgam_dtheta = ' symmat_print(model.dgam_dtheta) ';'];
    fprintf(f, '%s\n', str);
    
    %***************************
    %deta_dtheta (exogenous shocks)
    %***************************
    str = ['deta_dtheta = ' symmat_print(model.deta) ';'];
    fprintf(f, '%s\n', str);
    
    %***************************
    %deta_R (measurement error)
    %***************************
    str = ['dR_dtheta = ' symmat_print(model.dR) ';'];
    fprintf(f, '%s\n', str);
    fprintf(f, '%s\n', '');
else
    %Do not print the derivatives...and do not use those output arguments
end



fclose(f);

pause(.05); %A brief pause to ensure file has written to HD before accessing again.

%*******************************************************
% SYMMAT_PRINT:
% Print a symbolic expresison to a evaluateable function.
%*******************************************************
function str = symmat_print(x)
if isa(x, 'double')
    %Not actually symbolic
    str = mat2str(x);
elseif length(x(:)) == 1
    str = char(x);
else
    %Actually symbolic
    str = char(x);
    str = str(8:end-1);
    
    %Make into matlab matrix notation
    row_idx = findstr(str, '],');
    for j = 1:length(row_idx)
        str(row_idx(j):row_idx(j)+1)='];';
    end
    
end