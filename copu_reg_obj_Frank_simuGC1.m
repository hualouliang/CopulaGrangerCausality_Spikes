function [lk der]=copu_reg_obj_Frank_simuGC1(Y1,Y2,x,porder)

% copula GLM reduced model 1

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

mu1=x(1)+yhist(:,1:porder)*x(2:porder+1)';
mu2=x(porder+2)+yhist*x(porder+3:v-1)';

p1=exp(mu1)./(1+exp(mu1));
p2=exp(mu2)./(1+exp(mu2));
u1=1-p1;
u2=1-p2;

% theta=exp(x(end));
theta=x(end);

a1=frank_inv(u1,theta);
a2=frank_inv(u2,theta);
cr=frank_func(a1+a2,theta);
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

a1=frank_inv(u1,theta);
a2=frank_inv(u2,theta);
a4=frank_deri_t(a1+a2,theta);
a5=frank_func(a1+a2,theta);
a6=frankinv_deri_t(u1,theta); %% 
a7=frankinv_deri_t(u2,theta); %% 
a8=frank_deri_theta(a1+a2,theta);
a9=frankinv_deri_theta(u1, theta);
a10=frankinv_deri_theta(u2, theta);

%% Y1 =0 & Y2=0
indx10=(y1(:,1)==0);
indx20=(y2(:,1)==0);
indx00=indx10.*indx20;
indxtmp=find(indx00==1);

       fvec(indxtmp,1)=-a4(indxtmp)./a5(indxtmp).*p1(indxtmp).*(1-p1(indxtmp)).*a6(indxtmp);
       for i=2:porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,porder+2)=-a4(indxtmp)./a5(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)).*a7(indxtmp);
       for i=porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,porder+2).*yhist(indxtmp,i-porder-2);
       end       
%        fvec(indxtmp,v)=a8(indxtmp)./a5(indxtmp).*(a9(indxtmp)+a10(indxtmp))*exp(x(v));
       fvec(indxtmp,v)=a8(indxtmp)./a5(indxtmp).*(a9(indxtmp)+a10(indxtmp));
       
%% Y1 =0 & Y2=1       
indx10=(y1(:,1)==0);
indx21=(y2(:,1)==1);
indx01=indx10.*indx21;
indxtmp=find(indx01==1);

       fvec(indxtmp,1)=(a4(indxtmp).*p1(indxtmp).*(1-p1(indxtmp)).*a6(indxtmp)-p1(indxtmp).*(1-p1(indxtmp)))./(1-p1(indxtmp)-a5(indxtmp));
       for i=2:porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,porder+2)=(a4(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)).*a7(indxtmp))./(1-p1(indxtmp)-a5(indxtmp));
       for i=porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,porder+2).*yhist(indxtmp,i-porder-2);
       end         
%        fvec(indxtmp,v)=-a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))*exp(x(v))./(1-p1(indxtmp)-a5(indxtmp));
       fvec(indxtmp,v)=-a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))./(1-p1(indxtmp)-a5(indxtmp));
       

%% Y1 =1 & Y2=0       
indx11=(y1(:,1)==1);
indx20=(y2(:,1)==0);
indx10=indx11.*indx20;
indxtmp=find(indx10==1);

       fvec(indxtmp,1)=a4(indxtmp).*p1(indxtmp).*(1-p1(indxtmp)).*a6(indxtmp)./(1-p2(indxtmp)-a5(indxtmp));     
       for i=2:porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,porder+2)=(a4(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)).*a7(indxtmp)-p2(indxtmp).*(1-p2(indxtmp)))./(1-p2(indxtmp)-a5(indxtmp));  
       for i=porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,porder+2).*yhist(indxtmp,i-porder-2);
       end       
%        fvec(indxtmp,v)=-a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))*exp(x(v))./(1-p2(indxtmp)-a5(indxtmp));
       fvec(indxtmp,v)=-a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))./(1-p2(indxtmp)-a5(indxtmp));
       

%% Y1 =1 & Y2=1       
indx11=(y1(:,1)==1);
indx21=(y2(:,1)==1);
indx11=indx11.*indx21;
indxtmp=find(indx11==1);      

       fvec(indxtmp,1)=(p1(indxtmp).*(1-p1(indxtmp))-a4(indxtmp).*p1(indxtmp).*(1-p1(indxtmp)).*a6(indxtmp))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1);    
       for i=2:porder+1
           fvec(indxtmp,i)=fvec(indxtmp,1).*yhist(indxtmp,i-1);
       end
       fvec(indxtmp,porder+2)=(p2(indxtmp).*(1-p2(indxtmp))-a4(indxtmp).*p2(indxtmp).*(1-p2(indxtmp)).*a7(indxtmp))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1); 
       for i=porder+3:v-1
           fvec(indxtmp,i)=fvec(indxtmp,porder+2).*yhist(indxtmp,i-porder-2);
       end   
%        fvec(indxtmp,v)=a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))*exp(x(v))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1); 
       fvec(indxtmp,v)=a8(indxtmp).*(a9(indxtmp)+a10(indxtmp))./(p1(indxtmp)+p2(indxtmp)+a5(indxtmp)-1); 
       

%%       
fvec_new=fvec;

if nargout > 1
    der=-sum(fvec_new,1); 
end



end