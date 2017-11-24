function matrix2latex_black(matrix, filename, varargin)

% PURPOSE:  Export of a matlab array object for LaTeX, so that it can be
%           icluded as a tabular object with \input{out.tex}.
% ------------------------------------------------------------------------
% USAGE:    matrix2late(matrix, filename, varargs)
%
% WHERE:    matrix      = 2 dimensional numerical or cell array
%           filename    = valid filename, in which the resulting latex 
%                         code will be stored
%           varargs     = one ore more of the following (denominator, 
%                         value) combinations
%                       + 'rowLabels', array -> Can be used to label the 
%                       rows of the resulting latex table
%                       + 'columnLabels', array -> Can be used to label 
%                       the columns of the resulting latex table
%                           '...\\' = tile without any further cell content
%                           '\\'    = like \skip or \\ in text area
%                           '$..$'  = any math commands
%                       + 'alignment', 'value' -> Can be used to specify 
%                       the alginment of the table within the latex 
%                       document. Valid arguments are: 'l', 'c', and 'r' 
%                       for left, center, and right, respectively
%                       + 'format', 'value' -> Can be used to format the 
%                       input data. 'value' has to be a valid format 
%                       string, similar to the ones used in 
%                       fprintf('format', value);
%                       + 'size', 'value' -> One of latex' recognized 
%                       font-sizes, e.g. tiny, HUGE, Large, large, LARGE..
%
% NOTES:    Base programm was written by:
%           Author:   M. Koehler
%           Contact:  koehler@in.tum.de
%           Version:  1.1
%           Date:     May 09, 2004
% ------------------------------------------------------------------------
% RETURNS:  nothing since syntax is correct
% ------------------------------------------------------------------------
% EXAMPLE:
% matrix          = wfs;
% rowLabels       = {'Wald-F (lik)' 'P-Value' 'Wald-F (core)' 'P-Value'};
% columnLabels    = {'ltGap' 'hpGap' 'mvhpGap' 'stpGap' 'bpfoptGap' ...
%         'HKFGap' 'pfGap'};
% matrix2latex_table(matrix, 'U:\SNB\OUTPUT GAP\LaTeX\waldf.tex', ...
%     'rowLabels', rowLabels, 'columnLabels', columnLabels, 'alignment', ...
%     'c', 'format', '-6.4f', 'size', 'footnotesize');
%
% NOTES:    for row- and column-labels can also be used 
%           math-type-expressions like $\pi_{0}$.
% ------------------------------------------------------------------------
% SEE ALSO: matrix2latex_table, latex
% ========================================================================

% written by:
% Yves Blattmann
% Swiss National Bank
% 1. Departement / Economic Analysis
% Fraumünsterstrasse 8 
% CH-8022 Zürich
% yves.blattmann@snb.ch

    rowLabels = [];
    colLabels = [];
    alignment = 'l';
    format = [];
    fcst = [];
    textsize = [];
    if (rem(nargin,2) == 1 || nargin < 2)
        error('matrix2latex: ', 'Incorrect number of arguments to %s.', ...
            mfilename);
    end

    okargs = {'rowlabels','columnlabels', 'alignment', 'format', 'size', 'fcst'};
    for j=1:2:(nargin-2)
        pname = varargin{j};
        pval = varargin{j+1};
        k = strmatch(lower(pname), okargs);
        if isempty(k)
            error('matrix2latex: ', 'Unknown parameter name: %s.', pname);
        elseif length(k)>1
            error('matrix2latex: ', 'Ambiguous parameter name: %s.', ...
                pname);
        else
            switch(k)
                case 1  % rowlabels
                    rowLabels = pval;
                    if isnumeric(rowLabels)
                        rowLabels = cellstr(num2str(rowLabels(:)));
                    end
                case 2  % column labels
                    colLabels = pval;
                    if isnumeric(colLabels)
                        colLabels = cellstr(num2str(colLabels(:)));
                    end
                case 3  % alignment
                    alignment = lower(pval);
                    if alignment == 'right'
                        alignment = 'r';
                    end
                    if alignment == 'left'
                        alignment = 'l';
                    end
                    if alignment == 'center'
                        alignment = 'c';
                    end
                    if alignment ~= 'l' && alignment ~= 'c' && ...
                            alignment ~= 'r'
                        alignment = 'l';
                        warning('matrix2latex: ', ...
                            'Unkown alignment. (Set it to \''left\''.)');
                    end
                case 4  % format
                    format = lower(pval);
                case 5  % format
                    textsize = pval;
                case 6  % fcst
                    fcst = pval;
            end
        end
    end

    fid = fopen(filename, 'w');
    
    height = size(matrix, 1);
    if(~isempty(rowLabels))                 % rowLabels
        heightLabel = length(rowLabels);
    else                                    % no rowLabels
        heightLabel = height;
    end
    width = size(matrix, 2);

    if isnumeric(matrix)
        matrix = num2cell(matrix);
        for h=1:height
            for w=1:width
                if(~isempty(format))
                    matrix{h, w} = num2str(matrix{h, w}, format);
                else
                    matrix{h, w} = num2str(matrix{h, w});
                end
            end
        end
    end
    
    % center environment
    %fprintf(fid, '\\begin{center}\r\n');
    
    % ----- begin of style area ------------------------------------------
    
    % size
    if(~isempty(textsize))
        fprintf(fid, '\\begin{%s}\r\n', textsize);
    end
    
    % tabular environment
    fprintf(fid, '\t\\begin{tabular}{');
    
    if(~isempty(rowLabels))
        fprintf(fid, 'l');
    end
        
    for i=1:width
        fprintf(fid, '%c', alignment);
    end
    fprintf(fid, '}\r\n');
    fprintf(fid, '\t\\hline\r\n');
    
    % colon labels
    if(~isempty(colLabels))
        % the one of rowlabel is empty
        if(~isempty(rowLabels))
            fprintf(fid, '\t\t& ');
        end
        % others
        for w=1:width-1
            fprintf(fid, '%s & ', colLabels{w});  % \\textbf{%s}
        end
        % last
        fprintf(fid, '%s \\\\\r\n', ...           % \\textbf{%s}
            colLabels{width});
        % hline after colon labels
        fprintf(fid, '\t\t\\hline\r\n');
    end
    
    h = 1;
    for hL=1:heightLabel
        % write row labels
        if(~isempty(rowLabels))
            if (strcmp(rowLabels{hL},'\\'))
                 fprintf(fid, '\t\t\\\\\r\n', rowLabels{hL});                           
            elseif (~strcmp(rowLabels{hL}(end-1:end),'\\'))
                fprintf(fid, '\t\t%s & ', rowLabels{hL});            
                % write values
                for w=1:width-1
                    if(w<=width-fcst)
                        fprintf(fid, '%s & ', matrix{h, w}); % \\textbf{%s}
                    else
                        fprintf(fid, '%s & ', matrix{h, w});
                    end
                end
                fprintf(fid, '%s \\\\\r\n', ...
                    matrix{h, width});

                h = h + 1;
            else
                fprintf(fid, '\t\t\\textbf{%s} & ', rowLabels{hL}(1:end-2));
                % write empty &            
                for w=1:width-1
                   fprintf(fid, '& ');                  %, matrix{h, w});
                end
                fprintf(fid, '\\\\\r\n');               %, matrix{h, width});
            end
        end
    end
    
    % horizontal line at the end of tabel
    fprintf(fid, '\t\t\\hline\r\n');

    % ----- end of style area --------------------------------------------
    
    fprintf(fid, '\t\\end{tabular}\r\n');
    
    if(~isempty(textsize))
        fprintf(fid, '\\end{%s}\r\n', textsize);
    end
    
   % fprintf(fid,'\\end{center}\r\n');
    
    fclose(fid);
