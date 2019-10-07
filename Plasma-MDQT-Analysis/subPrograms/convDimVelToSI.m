function [velSI] = convDimVelToSI(velDim,hc,simParams)
%% Function notes
% This function converts velDim (a vector quantity in dimensionless units -
% e.g. in units of a*omega_pE) into a velSI (a vector quantity with units
% of m/s).

%% Function
% Calculate wigner seitz radius
a = (3/(4*pi*simParams.n))^(1/3);

% Calculate the ion's einstein frequency
w = sqrt(simParams.n*hc.e^2/(3*hc.mi*hc.eps));

% Convert velocity units
velSI = velDim.*a.*w;

end

