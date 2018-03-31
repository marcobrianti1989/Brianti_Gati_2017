function make_index(Y)
ny = length(Y);
for j=1:ny; 
    vn = char(Y(j));
    assignin('caller', [lower(vn), '_idx'], j);
end