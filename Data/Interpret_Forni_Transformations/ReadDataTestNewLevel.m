
function [data, codeData,tsc] = ReadDataTestNewLevel(transformations)

[x,text] = xlsread('DataTestNew.xlsx','Quarterly','B10:EA213');
if transformations == 0
    tsc = xlsread('DataTestNew.xlsx','Quarterly','B3:EA3');
elseif transformations == 1
    tsc = xlsread('DataTestNew.xlsx','Quarterly','B4:EA4');
end

%%% Data transformations
if transformations == 1
    for i=1:size(x,2)
        if tsc(i)==1
            data(:,i)=x(1:end,i);
            codeData(i)= 0;
        elseif tsc(i) == 2
            data(:,i)= x(:,i);
            codeData(i) = 0;
        elseif tsc(i)== 4
            data(:,i)=log(x(1:end,i))*100;
            codeData(i) = 0;
        elseif tsc(i) == 5
            data(:,i)=log(x(:,i))*100;
            codeData(i) = 0;
        end
    end
else
    for i=1:size(x,2)
        if tsc(i)==1
            data(:,i)=x(2:end,i);
            codeData(i)= 0;
        elseif tsc(i) == 2
            data(:,i)= diff(x(:,i));
            codeData(i) = 1;
        elseif tsc(i)== 4
            data(:,i)=log(x(2:end,i))*100;
            codeData(i) = 0;
        elseif tsc(i) == 5
            data(:,i)=diff(log(x(:,i)))*4;
            codeData(i) = 1;
            elseif tsc(i) == 6
            data(:,i)=diff(log(x(:,i)))*4;
            codeData(i) = 1;
        end
    end
end