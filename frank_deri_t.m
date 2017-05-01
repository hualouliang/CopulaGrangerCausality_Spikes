function tt=frank_deri_t(t, theta)
% derivation of Frank function
% Meng Hu @ Liang's lab at Drexel University, 2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

tmp=(1-exp(-theta)).*exp(-t);
tt=-1./theta.*(1./(1-tmp)).*tmp;

end