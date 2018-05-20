function make_index(Y,where_available)
% Input anything for the second argument if you wish the indexes to be
% available in the base function, not the caller.
ny = length(Y);
for j=1:ny; 
    vn = char(Y(j));
    
    switch nargin
        case 2
            assignin('caller', [lower(vn), '_idx'], j);
        case 1
            assignin('base', [lower(vn), '_idx'], j);
    end
end