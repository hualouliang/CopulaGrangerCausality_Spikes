function [gc12 gc21 para Lik Lik1 Lik2]=CopuReg_GC_Frank_fminunc(Y1,Y2,porder,options)


% Granger causality by copula GLM model

% Input:
% Y1 %% binary
% Y2 %% binary
% porder: model order
% options: parameter for fminunc
% % 

% Output:
% gc12: Granger causality from 1 to 2
% gc21: Granger causality from 2 to 1
% para: estimated parameters of full model
% Lik whole model
% Lik1 eliminate Y2->Y1 causal direction (binary->cont)
% Lik2 eliminate Y1->Y2 causal direction (cont->binary)

% Meng Hu @ Liang's lab at Drexel University, 08/2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

%% whole modeling

gamma_init=corr(Y1',Y2','type','Kendall');
gamma_init=copulaparam('Frank',gamma_init);
para_init=[rand(1,2+4*porder) gamma_init]; %% 1-4 for mu1; 5-8 for mu2; 9 for gamma
% tic
% options = optimset('GradObj','on','Display','iter','TolFun',1e-2,'TolX',1e-2,'LargeScale','off','MaxIter',200);
[x,fval]=fminunc(@(para) copu_reg_obj_Frank_simuGC(Y1,Y2,para,porder),para_init,options);
% toc        

para=x;
Lik=-fval;            
%% eliminate a causal direction (Y2->Y1) 
% flag2=0;

N=length(Y1);
para_init1=para_init([1:1+porder 2+2*porder:end]);
% tic
% options = optimset('GradObj','on','Display','iter','TolFun',1e-2,'TolX',1e-2,'LargeScale','off','MaxIter',200);
[x,fval]=fminunc(@(para) copu_reg_obj_Frank_simuGC1(Y1,Y2,para,porder),para_init1,options);
% toc        

para1=x;
Lik1=-fval;      

 %% eliminate a causal direction (Y1->Y2) (cont->binary)
% flag3=0;

N=length(Y1);
para_init2=para_init([1:2*porder+2 3*porder+3:end]);
% tic
% options = optimset('GradObj','on','Display','iter','TolFun',1e-2,'TolX',1e-2,'LargeScale','off','MaxIter',200);
[x,fval]=fminunc(@(para) copu_reg_obj_Frank_simuGC2(Y1,Y2,para,porder),para_init2,options);
% toc        

para2=x;
Lik2=-fval;      
  

%%

% Y2->Y1
gc21=real(Lik)-real(Lik1);

% Y1->Y2
gc12=real(Lik)-real(Lik2);


end