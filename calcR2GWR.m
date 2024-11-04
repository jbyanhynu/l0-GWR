function [R2,adj_R2,list_betas,AICc,totalLCN] = calcR2GWR(px,py,x,y,bw,pentalyCoe)
%bw:bandwidth
[row_px,~]=size(x);
list_influ=[];%hat matrix
list_betas=[];
totalLCN=zeros(row_px,1);
parfor k2=1:row_px
    distance_matrix=pdist2([px,py],[px(k2),py(k2)],'euclidean');
    temp_distance_matrix=sort(distance_matrix);
    bandWidth=temp_distance_matrix(bw)*1.0000001;
    kernel=(1 - (distance_matrix/bandWidth).^2).^2;%kernel:bi-square
    kernel(distance_matrix>=bandWidth)=0;%kernel:bi-square
    %kernel=exp(-(distance_matrix/bandWidth));%exp kernel
%     kernel=exp(-1/2*((distance_matrix/bandWidth).^2));%Gauss
%     kernel(distance_matrix>=bandWidth)=0;
    Xw=x.*(sqrt(kernel));
    Yw=y.*(sqrt(kernel));
    betas = ABESS(Xw,Yw,pentalyCoe,bw);
    tempbeta=betas';
    position=(tempbeta~=0);
    xxx=Xw(:,position);
    %LCN
    xw=xxx;
    sxw=sqrt(sum(xw.^2));
    temp=(xw'./sxw');
    sxw=temp';
    svdx=svd(sxw);
    totalLCN(k2)=svdx(1)/svdx(end);
    
    xtx_inv_xt=(xxx'*xxx)\xxx';
    influ=xxx(k2,:)*xtx_inv_xt(:,k2);
    list_betas=[list_betas;betas];
    list_influ=[list_influ,influ];
end
predy=sum(list_betas.*x,2);
list_resid=y-predy;
%y_hat=S*y��
tr_S=sum(list_influ);
% %�в�ƽ����
SSE=list_resid'*list_resid;
% %���ƽ����
SST=(y-mean(y))'*(y-mean(y));
% %�������Ŷ��Լ�����֮�������Ŷ�
R2=1-SSE/SST;
% %�������ɶȵ���֮���R2,�ų��Ա��������뵼�·�������Ŷ���ǿ
adj_R2=1 - (1 - R2) * (row_px - 1) / ((row_px) - tr_S - 1);
sigma_hat=sqrt(SSE/row_px);%normalised residual sum of squares from the local regression and is defined as
% matlab��log��ʾ��������Ȼ����Ϊ�׵Ķ���
AICc=2*row_px*log(sigma_hat)+row_px*log(2*pi)+row_px*((row_px+tr_S)/(row_px-2-tr_S));
end