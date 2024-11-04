function [y,px,py,x] = readreal(path)
data=xlsread(path,1);
y=data(:,1);
px=data(:,2);
py=data(:,3);
x=data(:,4:end);
intercept=ones(length(px),1);
x=(x-mean(x))./std(x);
x=[intercept,x];
y=(y-mean(y))./std(y);
end

