function tt=frank_func(t,theta)
% Frank function
% Meng Hu @ Liang's lab at Drexel University, 2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

tt=-1./theta.*log(1-(1-exp(-theta)).*exp(-t));

end