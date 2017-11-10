% IR_FIG - Generate "nice" impulse repsonse figures from a matrix of IR
% data
%
% usage
%
% fig_out = ir_fig(y,var_list, fig_title,lstyle(optional))
%
% where
%
% y = a TxN matrix of impulse responses
% var_list = name of the variables plotted
% fig_title = Title for entire figure
% lstyle = a line style (eg 'r-*' means red, stared line). '-' is default.



function fig = ir_fig(y, var_list, fig_title, varargin)

%Create figure, set size
fig = figure;
set(fig, 'PaperOrientation','landscape', 'PaperPosition', [.25, .25, 11,8]);


%Figure dimensions
n = length(var_list);
n_row = 2;
n_col = ceil(n/n_row);
n_slide = size(y,3);  %Are the confidence bands
if size(y,2) ~= n
    error('Unequal number of columns and variable titles');
end
t = size(y,1);

%Choose Line Style
if ~isempty(varargin)
    lstyle = varargin{1};
else
    lstyle = {'-'};
end

%ylim(1) = min(min(min(y))); ylim(1)= ylim(1) - .1*max(abs(ylim));
%ylim(2) = max(max(max(y))); ylim(2)= ylim(2) + .1*max(abs(ylim));
%Plot each IR
for j = 1:n
    s = subplot(n_row,n_col,j);
    hold on
    for k = 1:n_slide
        p = plot([1:t]', y(:,j,k), lstyle{k}, 'linewidth', 1);
    end 
    plot([1:t]', 0*[1:t]', ':k');
    title(var_list{j},'interpreter', 'latex', 'fontsize', 12);
    set(s, 'xlim', [1,t], 'ylim', ylim );
    if j ==1
        xlabel('period')
        ylabel('pct deviation from ss');
    end
end

%Title on figure
if ~isempty(fig_title)
    fig_title
   suptitle(fig,fig_title, 'tex', 16);
end

%Legend
if n_slide > 1 && ~isempty(varargin{2})
   
 
     l=legend(varargin{2}, 'Location', 'SouthEast');
    %set(l, 'OuterPosition', [0.63 0.15 0.219 0.248], 'FontSize', 14,'Interpreter', 'tex')

    %set(l, 'location', 'southeastoutside');
    %set(l,'interpreter', 'latex', 'fontsize', 10, 'string',varargin{2});
end

function hout=suptitle(fig, str, varargin)
%SUPTITLE Puts a title above all subplots.
%	SUPTITLE('text', interpreter (optional), fontsize (optional)) 
%   adds text to the top of the figure
%	above all subplots (a "super title"). Use this function
%	after all subplot commands.

% Drea Thomas 6/15/95 drea@mathworks.com
% Modified by Ryan Chahrour 4-16-2008 to allow interpreter option

% Warning: If the figure or axis units are non-default, this
% will break.

% Parameters used to position the supertitle.

figure(fig);
% Amount of the figure window devoted to subplots
plotregion = .92;

% Y position of title in normalized coordinates
titleypos  = .95;

% Fontsize for supertitle
fs = get(gcf,'defaultaxesfontsize')+4;

% Fudge factor to adjust y spacing between subplots
fudge=1;

haold = gca;
figunits = get(gcf,'units');

% Get the (approximate) difference between full height (plot + title
% + xlabel) and bounding rectangle.

	if (~strcmp(figunits,'pixels')),
		set(gcf,'units','pixels');
		pos = get(gcf,'position');
		set(gcf,'units',figunits);
    else
		pos = get(gcf,'position');
	end
	ff = (fs-4)*1.27*5/pos(4)*fudge;

       % The 5 here reflects about 3 characters of height below
       % an axis and 2 above. 1.27 is pixels per point.

% Determine the bounding rectange for all the plots

h = findobj(fig,'Type','axes');
max_y=0;
min_y=1;

oldtitle =0;
for i=1:length(h),
	if (~strcmp(get(h(i),'Tag'),'suptitle')),
		pos=get(h(i),'pos');
		if (pos(2) < min_y), min_y=pos(2)-ff/5*3;end;
		if (pos(4)+pos(2) > max_y), max_y=pos(4)+pos(2)+ff/5*2;end;
	else,
		oldtitle = h(i);
	end
end

if max_y > plotregion,
	scale = (plotregion-min_y)/(max_y-min_y);
	for i=1:length(h),
		pos = get(h(i),'position');
		pos(2) = (pos(2)-min_y)*scale+min_y;
		pos(4) = max(.01,pos(4)*scale-(1-scale)*ff/5*3);
		set(h(i),'position',pos);
	end
end

np = get(gcf,'nextplot');
set(gcf,'nextplot','add');
if (oldtitle),
	%delete(oldtitle);
end
ha=axes('pos',[0 1 1 1],'visible','off','Tag','suptitle');
ht=text(.5,titleypos-1,str);
set(ht,'horizontalalignment','center','fontsize',fs);
set(gcf,'nextplot',np);
axes(haold);
if nargout,
	hout=ht;
end

if ~isempty(varargin)
set(ht, 'interpreter', varargin{1}, 'fontsize', varargin{2});
end
