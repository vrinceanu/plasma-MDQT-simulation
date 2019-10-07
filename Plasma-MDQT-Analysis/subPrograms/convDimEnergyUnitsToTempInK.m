function [matrixInK] = convDimEnergyUnitsToTempInK(matrixInDimUnits,hc,simParams)
%% Function Notes
% This function takes a matrix input that has energy in units of Ec =
% e^2/(4*pi*eps*a)) and converts it to an energy in units of temperature (K).

%% Calculate Wigner Sitz Radius
a = (3/(4*pi*simParams.n))^(1/3);

%% Convert from dim units to K
matrixInK = matrixInDimUnits.*hc.e^2/(4*pi*hc.eps*hc.kb*a);

end

