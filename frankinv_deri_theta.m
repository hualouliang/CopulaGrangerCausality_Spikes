function tt=frankinv_deri_theta(t, theta)
% inverse Frank function's derivation according to theta
% Meng Hu @ Liang's lab at Drexel University, 2014
% Paper: "Copula Regression Analysis of Simultaneously Recorded 
% Frontal Eye Field and Inferotemporal Spiking Activity During Object-based Working Memory"

tt=t.*exp(-theta*t)./(exp(-theta*t)-1)-exp(-theta)./(exp(-theta)-1);

end