function tt=frank_deri_theta(t, theta)
% Frank function derivation according to theta
% Meng Hu @ Liang's lab at Drexel University, 2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

tmp=1-(1-exp(-theta)).*exp(-t);
tt=1./theta^2.*log(tmp)+1./theta./tmp.*exp(-t).*exp(-theta);

end