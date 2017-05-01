function [lk der]=copu_reg_obj_Gauss_simuGC2(Y1,Y2,x,porder)

% Gaussian copula GLM reduced model 2

% Input:
% Y1 %% binary
% Y2 %% binary
% porder: model order
% x: initiative parameter
% % 

% Output:
% lk %% log-likelihood
% der %% derivation

% Meng Hu @ Liang's lab at Drexel University, 08/2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"


N=length(Y1);
v=length(x);

win=1;
y1(:,1)=Y1(1+porder*win:end);
y2(:,1)=Y2(1+porder*win:end);
ll=length(y1);

for n=1:porder
    for m=1:ll
        y1(m,n+1)=sum(Y1(porder*win+m-1-(n-1)*win:-1:porder*win+m-n*win));
        y2(m,n+1)=sum(Y2(porder*win+m-1-(n-1)*win:-1:porder*win+m-n*win));
    end
end
yhist=[y1(:,2:end) y2(:,2:end)];

mu1=x(1)+yhist*x(2:2*porder+1)';
mu2=x(2*porder+2)+yhist(:,porder+1:end)*x(2*porder+3:v-1)';

p1=exp(mu1)./(1+exp(mu1));
p2=exp(mu2)./(1+exp(mu2));
u1=1-p1;
u2=1-p2;

gamma=-1+2*exp(-exp(x(v)));
sigma=[1 gamma ; gamma 1];

a1=norminv(u1,0,1);
a2=norminv(u2,0,1);
% a4=mvnpdf([a1 a2],[0 0],sigma); 
cr=mvncdf([a1 a2],[0 0],sigma);
LLK=[];

%% Y1 =0 & Y2=0
indx10=(y1(:,1)==0);
indx20=(y2(:,1)==0);
indx00=indx10.*indx20;
indxtmp=find(indx00==1);

LLK(indxtmp)=cr(indxtmp);
%% Y1 =0 & Y2=1       
indx10=(y1(:,1)==0);
indx21=(y2(:,1)==1);
indx01=indx10.*indx21;
indxtmp=find(indx01==1);

LLK(indxtmp)=1-p1(indxtmp)-cr(indxtmp);
%% Y1 =1 & Y2=0       
indx11=(y1(:,1)==1);
indx20=(y2(:,1)==0);
indx10=indx11.*indx20;
indxtmp=find(indx10==1);

LLK(indxtmp)=1-p2(indxtmp)-cr(indxtmp);
%% Y1 =1 & Y2=1       
indx11=(y1(:,1)==1);
indx21=(y2(:,1)==1);
indx11=indx11.*indx21;
indxtmp=find(indx11==1);      

LLK(indxtmp)=p1(indxtmp)+p2(indxtmp)+cr(indxtmp)-1;

%%
lk=sum(log(LLK));
lk=-lk;


%% derivation

gamma=-1+2*exp(-exp(x(v)));
sigma=[1 gamma ; gamma 1];

a1=norminv(u1,0,1);
a2=norminv(u2,0,1);
% a4=mvnpdf([a1 a2],[0 0],sigma); 
a5=mvncdf([a1 a2],[0 0],sigma);
a6=normcdf((a2-gamma*a1)./sqrt(1-gamma^2),0,1);
a7=normcdf((a1-gamma*a2)./sqrt(1-gamma^2),0,1);

%% Y1 =0 & Y2=0
indx10=(y1(:,1)==0);
indx20=(y2(:,1)==0);
indx00=indx10.*indx20;
indxtmp=find(indx00==1);

       fvec(indxtmp,1)=-a6(indxtmp)./a5(indxtmp).*p1(indxtmp).*(1-p1(indxtmp));
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,2*porder+2)=-a7(indxtmp)./a5(indxtmp).*p2(indxtmp).*(1-p2(indxtmp));
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-porder-2);
       end       
       fvec(indxtmp,v)=1./a5(indxtmp).*exp(-(a1(indxtmp).^2-2*a1(indxtmp).*a2(indxtmp)*gamma+a2(indxtmp).^2)/(2*(1-gamma.^2)))/(2*pi*sqrt(1-gamma.^2))*(-2*exp(-exp(x(v)))*exp(x(v)));

%% Y1 =0 & Y2=1       
indx10=(y1(:,1)==0);
indx21=(y2(:,1)==1);
indx01=indx10.*indx21;
indxtmp=find(indx01==1);

       fvec(indxtmp,1)=(a6(indxtmp).*p1(indxtmp).*(1-p1(indxtmp))-p1(indxtmp).*(1-p1(indxtmp)))./(1-p1(indxtmp)-a5(indxtmp));
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,2*porder+2)=(a7(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)))./(1-p1(indxtmp)-a5(indxtmp));
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-porder-2);
       end       
       fvec(indxtmp,v)=-1./(1-p1(indxtmp)-a5(indxtmp)).*exp(-(a1(indxtmp).^2-2*a1(indxtmp).*a2(indxtmp)*gamma+a2(indxtmp).^2)/(2*(1-gamma.^2)))/(2*pi*sqrt(1-gamma.^2))*(-2*exp(-exp(x(v)))*exp(x(v)));

%% Y1 =1 & Y2=0       
indx11=(y1(:,1)==1);
indx20=(y2(:,1)==0);
indx10=indx11.*indx20;
indxtmp=find(indx10==1);

       fvec(indxtmp,1)=a6(indxtmp).*p1(indxtmp).*(1-p1(indxtmp))./(1-p2(indxtmp)-a5(indxtmp));     
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end 
       fvec(indxtmp,2*porder+2)=(a7(indxtmp).*p2(indxtmp).*(1-p2(indxtmp))-p2(indxtmp).*(1-p2(indxtmp)))./(1-p2(indxtmp)-a5(indxtmp)); 
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-porder-2);
       end      
       fvec(indxtmp,v)=-1./(1-p2(indxtmp)-a5(indxtmp)).*exp(-(a1(indxtmp).^2-2*a1(indxtmp).*a2(indxtmp)*gamma+a2(indxtmp).^2)/(2*(1-gamma.^2)))/(2*pi*sqrt(1-gamma.^2))*(-2*exp(-exp(x(v)))*exp(x(v)));

%% Y1 =1 & Y2=1       
indx11=(y1(:,1)==1);
indx21=(y2(:,1)==1);
indx11=indx11.*indx21;
indxtmp=find(indx11==1);      

       fvec(indxtmp,1)=(p1(indxtmp).*(1-p1(indxtmp))-a6(indxtmp).*p1(indxtmp).*(1-p1(indxtmp)))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1);  
       for i=2:2*porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,2*porder+2)=(p2(indxtmp).*(1-p2(indxtmp))-a7(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1); 
       for i=2*porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,2*porder+2).*yhist(indxtmp,i-porder-2);
       end   
       fvec(indxtmp,v)=1./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1).*exp(-(a1(indxtmp).^2-2*a1(indxtmp).*a2(indxtmp)*gamma+a2(indxtmp).^2)/(2*(1-gamma.^2)))/(2*pi*sqrt(1-gamma.^2))*(-2*exp(-exp(x(v)))*exp(x(v)));

%%       
fvec_new=fvec;

if nargout > 1
    der=-sum(fvec_new,1); 
end


end