function Index=ExtractIndices(Index)
Index=Index(:);
T=length(Index);
index=[];
for tt=1:T
    if Index(tt)~=0
        index=[index; tt];
    end
end
Index=index;




