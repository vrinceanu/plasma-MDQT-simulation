function [temps] = solveRateEquations(c,times,k,T0)
%% Program notes

%   This program solves the system of differential equations described in OneNote entry 04/11/19-04/12/19: Laser Cooling 
%   Simulations - Correlation Temperature. This is a system of equations that describes plasma temperature evolution
%   during 1D optical molasses. The system of equations solves for the temperature parallel/perpendicular to the laser
%   cooling axis and the correlation energy.

%   Inputs
%       c: (1x3 double) rate coefficients to be used in the model [beta nu mu]. beta is the laser cooling rate. nu is
%                       the cross-axis thermalization rate. mu is the rate that the correlation energy relaxes to 
%                       equilibrium. The units for each quantity should be us^-1
%       times: (nx1 double) column vector containing time points for which the temperatures should be solved at in units
%                           of us
%       k: (1x1 double) plasma screening parameter
%       T0: (1x3 double) initial plasma temperatures [T_|| T_perp Uii/kb] in dimensionless units

%   Outputs
%       (3nx1 double) column vector containing the solutions for [T_||; T_perp; Uii/kb] in dimensionless units
%                   *note* that this is 3x as long as 'times'. The temperatures are intentionally 'stacked' this way so
%                   that the output of this function is compatible with 'nlinfit'.

%% Solve the rate equations
% Define the system of rate equations to be solved
odefun = @(t,T) defineRateEqns(t,T,c,k);
% Solve the rate equations
[times,temps] = ode45(odefun,times,T0);
% Stack temperature solutions so that function output is compatible with 'nlinfit'
temps = temps(:);

end

