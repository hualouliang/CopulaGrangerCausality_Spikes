function tt=frank_inv(t,theta)
% inverse Frank function
% Meng Hu @ Liang's lab at Drexel University, 2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

tt=-log((exp(-theta.*t)-1)./(exp(-theta)-1));

end