% The data is generated 
% for comparison glm and Gaussian Copula

function [dat] = gendata_gc(M,N,rho)

% M: number of trial
% N: number of sample points
% rho: copula coefficient
% beta: parameters for marignal regression

% Meng Hu @ Liang's lab at Drexel University, 08/2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"


beta1 = 0.5; beta2 = -0.5; beta3=0;
beta4 = -2; beta5 = 1; beta6=0.5; %% Y1->Y2

sigma = repmat(rho,2,2);
sigma = sigma-diag(diag(sigma))+eye(2);

dat=[];
for i=1:M    
    
    Y1=[];
    Y2=[];
    Y1(1)=1;
    Y2(1)=1;
    U = copularnd('Gaussian',sigma,N); %    
    for j=2:N
        p1 = 1 ./ (1+ exp(-(beta1 + beta2*Y1(j-1) + beta3*Y2(j-1)))); % x2 ~ N(0, 1)
        p2 = 1 ./ (1+ exp(-(beta4 + beta5*Y1(j-1) + beta6*Y2(j-1)))); % x2 ~ N(0, 1)
        Y1(j) = U(j,1) < p1; 
        Y2(j) = U(j,2) < p2;    
    end
    dat(i,:,:) = [Y1' Y2'];
    
end


end
