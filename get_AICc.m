function AICc = get_AICc(bw,px,py,x,y,pentalyCoe)% aims to get the minimum AICc
row_px=size(x,1);
bw=round(bw);
list_influ=zeros(1,row_px);%freedom of matrix
list_betas=[];
parfor k2=1:row_px
    distance_matrix=pdist2([px,py],[px(k2),py(k2)],'euclidean');% calculating the weights matrix
    temp_distance_matrix=sort(distance_matrix);% sorting of the distance
    bandWidth=temp_distance_matrix(bw)*1.0000001;
    %% ÆäÓàºË
    kernel=(1 - (distance_matrix/bandWidth).^2).^2;%kernel:bi-square
    kernel(distance_matrix>=bandWidth)=0;%kernel:bi-square
    %kernel=exp(-(distance_matrix/bandWidth));%exp kernel
	%kernel=exp(-1/2*((distance_matrix/bandWidth).^2));%gauss kernel
	%kernel(distance_matrix>=bandWidth)=0;
    %cite from this paper: 
    %<simultaneous coefficient penalization and model selection in geographically weighted regression : the geographically weighted lasso>
    Xnew=x.*(sqrt(kernel));
    Ynew=y.*(sqrt(kernel));
    betas = ABESS(Xnew,Ynew,pentalyCoe,bw);
    tempbeta=betas';
    position=(tempbeta~=0);
    xxx=Xnew(:,position);
    xtx_inv_xt=(xxx'*xxx)\xxx';
    influ=xxx(k2,:)*xtx_inv_xt(:,k2);
    list_betas=[list_betas;betas];
    list_influ(k2)=influ;
end
predy=sum(list_betas.*x,2);
list_resid=y-predy;
tr_S=sum(list_influ);
SSE=list_resid'*list_resid;
sigma_hat=sqrt(SSE/row_px);
AICc=2*row_px*log(sigma_hat)+row_px*log(2*pi)+row_px*((row_px+tr_S)/(row_px-2-tr_S));
end