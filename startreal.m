clear;
tic;
%��1��������
[y,px,py,x]=readreal('./GData_utm.csv');
[~,p]=size(x);
qmin=40+2*p;
qmax=length(px);

%��2����ȡ�ͷ����������������ο���ѡ��1.�Ż��ķ�ʽ��pentalyMin��pentalyMax����2.�̶�ֵ��pentalyCoe=log(log(p));
%1. �Ż��ķ�ʽ��(0.5log(log(p)),2*log(p)��;��log(log(p)),log(p)��
%options = optimset('fminbnd');
%options.Display='iter';
%pentalyMin=log(log(p));%raw:0.5
%pentalyMax=log(p);%raw:2
%pentalyCoe=fminbnd('pentalyL',pentalyMin,pentalyMax,options,px,py,x,y);%ͨ��˫ѭ����ȡ�����ֵ

%2.�̶�ֵ�ķ�ʽ
 pentalyCoe=log(log(p));

%��3������pentalyCoe��ȡ����
options = optimset('fminbnd');
options.Display='iter';
band_width=fminbnd('get_AICc',qmin,qmax,options,px,py,x,y,pentalyCoe);
band_width=round(band_width);
%��4����������
[R2,adjR2,Finalbetas,AICc,totalLCN]=calcR2GWR(px,py,x,y,band_width,pentalyCoe);
tt=Finalbetas;
[numberSample,colx]=size(tt);
positiveNumber=sum(tt>0)/numberSample*100;%��Ԫ�صĸ���
zeroNumber=sum(tt==0)/numberSample*100;%��Ԫ�صĸ���
negativeNumber=sum(tt<0)/numberSample*100;%��Ԫ�صĸ���
ics=sum(abs(positiveNumber-negativeNumber)+zeroNumber)/colx;%����ģ�͵�ϵ�������׽�����