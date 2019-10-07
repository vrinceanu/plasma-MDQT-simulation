function [dTdt] = defineRateEqns(t,T,c,k)
%% Program notes

%   This program defines the system of differential equations described in OneNote entry 04/11/19-04/12/19: Laser Cooling 
%   Simulations - Correlation Temperature. This is a system of equations that describes plasma temperature evolution
%   during 1D optical molasses. The system of equations solves for the temperature parallel/perpendicular to the laser
%   cooling axis and the correlation energy.

%   Inputs
%       t: (1x1 double) 'ode45' requires the function handle to accept a time input, but it goes unused
%       T: (1x3 double) plasma temperature in dimensionless units [T_|| T_perp Uii/kb]
%       c: (1x3 double) rate equation coefficients with units us^-1 [beta nu mu]
%       k: (1x1 double) plasma screening parameter

%   Outputs
%       dTdt: (3x1 double) contains the instantaneous rate of change of T in dimensionless units

%% Define differential equations
% Calculate average plasma temperature/coupling parameter for use in correlation energy calculation
Tav = (T(1)+2*T(2))/3;
gam = 1/Tav;
% Calculate correlation energy 
UiiEQ = calcUiiEQ(gam,k);
% Define rate equations
dTdt = zeros(3,1);
dTdt(3) = c(3)*(UiiEQ-T(3));
dTdt(2) = -c(2)*(T(2)-T(1))-2/3*dTdt(3);
dTdt(1) = -2*c(1)*T(1)+2*c(2)*(T(2)-T(1))-2/3*dTdt(3);

end

