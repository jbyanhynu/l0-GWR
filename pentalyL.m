function [AICc,band_width,R2,adj_R2,list_betas] = pentalyL(pentalyCoe,px,py,x,y)
[~,p]=size(x);
qmin=40+2*p;
qmax=length(px);
options2 = optimset('fminbnd');
options2.TolX=1e0;
band_width=fminbnd('get_AICc',qmin,qmax,options2,px,py,x,y,pentalyCoe);%获得当前岭值下的最佳带宽
[R2,adj_R2,list_betas,AICc]=calcR2GWR(px,py,x,y,round(band_width),pentalyCoe);
end