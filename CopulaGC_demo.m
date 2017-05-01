% Demo codes for the Copula-based Granger causality for spike train data
% 
% Meng Hu @ Liang's lab at Drexel University, 08/2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"


clear
%% Simulated data generation

% number of trials for the simulated data
Ntrial=50;

% Length of simulated data
Npoint=1000;

% copula coefficient for generating simulated data
rho=0.5;


% generate bivariate data with causal influence
% marginal regressions:
% Y1: logit(p1) = beta1 + beta2*Y1(t-1) + beta3*Y2(t-1)
% Y2: logit(p2) = beta4 + beta5*Y1(t-1) + beta6*Y2(t-1)
% 
% beta1 = 0.5
% beta2 = -0.5
% beta3 = 0
% beta4 = -2
% beta5 = 1
% beta6 = 0.5
% 
% Thus, true Granger causality is Y1 -> Y2

dat = gendata_gc(Ntrial,Npoint,rho);

% model order
porder=1;

% parameter for model estimation
options = optimset('GradObj','on','Display','off','TolFun',1e-4,'TolX',1e-4,'LargeScale','off','MaxIter',200);


%% Copula Granger causality (Frank)

gc12_frank=[];
gc21_frank=[];
para_frank=[];
for n=1:size(dat,1)
    Y1=squeeze(dat(n,:,1));
    Y2=squeeze(dat(n,:,2));
% designed Granger causality: Y1->Y2
    try
% Frank copula
        [gc12_frank(n) gc21_frank(n) para_frank(:,n)]=CopuReg_GC_Frank_fminunc(Y1,Y2,porder,options);
    end  
end

%%%%%%%%% Comparing the true and estimated parameters

beta_est = mean(para_frank(1:6,:)');
ken_est = [];
for i = 1:Ntrial
%     rho_est=[rho_est,copulaparam('Gaussian',copulastat('Frank',para_frank(7,i)))];
    ken_est=[ken_est,copulastat('Frank',para_frank(7,i))];
end
ken_est = mean(ken_est);

para_true = [0.5 -0.5 0 -2 1 0.5 copulastat('Gaussian',0.5)];
para_est = [beta_est ken_est];

figure
bar([para_true; para_est]')
xlabel('Parameter estimation')
set(gca,'XTick',[1:7],'XTicklabel',{' Beta1 ' ' Beta2 ' ' Beta3 ' ' Beta4 ' ' Beta5 ' ' Beta6 ' 'Kendal'});
ax=axis;
axis([ax(1) ax(2) -2.2 1.2])
legend('True','Estimated')


%%%%%%%%% Granger causality

dat_perm1 = cat(3,dat(randperm(Ntrial),:,1),dat(randperm(Ntrial),:,2));
dat_perm2 = cat(3,dat(randperm(Ntrial),:,1),dat(randperm(Ntrial),:,2));
dat_perm = cat(1,dat_perm1,dat_perm2);

gc12_frank_perm=[];
gc21_frank_perm=[];
for n=1:size(dat_perm,1)
    Y1_perm=squeeze(dat_perm(n,:,1));
    Y2_perm=squeeze(dat_perm(n,:,2));
% designed Granger causality: Y1->Y2
    try
% Frank copula
        [gc12_frank_perm(n) gc21_frank_perm(n)]=CopuReg_GC_Frank_fminunc(Y1_perm,Y2_perm,porder,options);
    end  
end

sig_alpha = 0.05;

gc12tmp = sort(gc12_frank_perm);
gc12_thresh = gc12tmp(fix(length(gc12tmp)*(1-sig_alpha)));
gc21tmp = sort(gc21_frank_perm);
gc21_thresh = gc21tmp(fix(length(gc21tmp)*(1-sig_alpha)));

figure;
bar([mean(gc12_frank) mean(gc21_frank)])
hold on
plot([1,2],[gc12_thresh,gc21_thresh],'color','red','marker','^')
set(gca,'XTick',[1:2],'XTicklabel',{'Y1 -> Y2' 'Y2 -> Y1'});
axis([0 3 0 18.5])
legend('Estimated GC','95% significance level')
ylabel('Granger causality')


%%  Copula Granger causality (Gaussian); Note the Gaussian code would run slower than the Frank


gc12_Gauss=[];
gc21_Gauss=[];
para_Gauss=[];
for n=1:size(dat,1)
    Y1=squeeze(dat(n,:,1));
    Y2=squeeze(dat(n,:,2));
% designed Granger causality: Y1->Y2
    try
% Frank copula
        [gc12_Gauss(n) gc21_Gauss(n) para_Gauss(:,n)]=CopuReg_GC_Gauss_fminunc(Y1,Y2,porder,options);
    end  
    
    n
end

%%%%%%%%% Comparing the true and estimated parameters

beta_est = mean(para_Gauss(1:6,:)');
ken_est = [];
for i = 1:Ntrial
%     rho_est=[rho_est,copulaparam('Gaussian',copulastat('Frank',para_frank(7,i)))];
    tmpgau = -1+2*exp(-exp(para_Gauss(7,i)));
    ken_est=[ken_est,copulastat('Gaussian',tmpgau)];
end
ken_est = mean(ken_est);

para_true = [0.5 -0.5 0 -2 1 0.5 copulastat('Gaussian',0.5)];
para_est = [beta_est ken_est];

figure
bar([para_true; para_est]')
xlabel('Parameter estimation')
set(gca,'XTick',[1:7],'XTicklabel',{' Beta1 ' ' Beta2 ' ' Beta3 ' ' Beta4 ' ' Beta5 ' ' Beta6 ' 'Kendal'});
ax=axis;
axis([ax(1) ax(2) -2.2 1.2])
legend('True','Estimated')


%%%%%%%%% Granger causality

dat_perm1 = cat(3,dat(randperm(Ntrial),:,1),dat(randperm(Ntrial),:,2));
dat_perm2 = cat(3,dat(randperm(Ntrial),:,1),dat(randperm(Ntrial),:,2));
dat_perm = cat(1,dat_perm1,dat_perm2);

gc12_Gauss_perm=[];
gc21_Gauss_perm=[];
for n=1:size(dat_perm,1)
    Y1_perm=squeeze(dat_perm(n,:,1));
    Y2_perm=squeeze(dat_perm(n,:,2));
% designed Granger causality: Y1->Y2
    try
% Frank copula
        [gc12_Gauss_perm(n) gc21_Gauss_perm(n)]=CopuReg_GC_Gauss_fminunc(Y1_perm,Y2_perm,porder,options);
    end  
end

sig_alpha = 0.05;

gc12tmp = sort(gc12_Gauss_perm);
gc12_thresh = gc12tmp(fix(length(gc12tmp)*(1-sig_alpha)));
gc21tmp = sort(gc21_Gauss_perm);
gc21_thresh = gc21tmp(fix(length(gc21tmp)*(1-sig_alpha)));

figure;
bar([mean(gc12_Gauss) mean(gc21_Gauss)])
hold on
plot([1,2],[gc12_thresh,gc21_thresh],'color','red','marker','^')
set(gca,'XTick',[1:2],'XTicklabel',{'Y1 -> Y2' 'Y2 -> Y1'});
axis([0 3 0 18.5])
legend('Estimated GC','95% significance level')
ylabel('Granger causality')



