clear;
tic;
%（1）读数据
[y,px,py,x]=readreal('./GData_utm.csv');
[~,p]=size(x);
qmin=40+2*p;
qmax=length(px);

%（2）获取惩罚参数。有两种情形可以选择：1.优化的方式【pentalyMin，pentalyMax】；2.固定值：pentalyCoe=log(log(p));
%1. 优化的方式【(0.5log(log(p)),2*log(p)】;【log(log(p)),log(p)】
%options = optimset('fminbnd');
%options.Display='iter';
%pentalyMin=log(log(p));%raw:0.5
%pentalyMax=log(p);%raw:2
%pentalyCoe=fminbnd('pentalyL',pentalyMin,pentalyMax,options,px,py,x,y);%通过双循环获取最佳岭值

%2.固定值的方式
 pentalyCoe=log(log(p));

%（3）依据pentalyCoe获取带宽
options = optimset('fminbnd');
options.Display='iter';
band_width=fminbnd('get_AICc',qmin,qmax,options,px,py,x,y,pentalyCoe);
band_width=round(band_width);
%（4）参数估计
[R2,adjR2,Finalbetas,AICc,totalLCN]=calcR2GWR(px,py,x,y,band_width,pentalyCoe);
tt=Finalbetas;
[numberSample,colx]=size(tt);
positiveNumber=sum(tt>0)/numberSample*100;%正元素的个数
zeroNumber=sum(tt==0)/numberSample*100;%零元素的个数
negativeNumber=sum(tt<0)/numberSample*100;%负元素的个数
ics=sum(abs(positiveNumber-negativeNumber)+zeroNumber)/colx;%整个模型的系数符号易解释性