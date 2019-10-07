function [UiiEQ] = calcUiiEQ(gam,k)
%% Program notes

%   This program calculates the UNP correlation energy as described in OneNote entry 05/06/19: Troubleshoot Issues with 
%   MDQT Code v9.

%   Inputs
%       gam: (1x1 double) average plasma coupling parameter
%       k: (1x1 double) plasma screening parameter

%   Outputs
%       UiiEQ: (1x1 double) plasma correlation energy in dimensionless units


%% Calculate correlation energy
% See the referenced OneNote entry for more information.
Ebcc = -0.895929-0.103731*k^2+0.003084*k^4-0.000131*k^6;
da = -0.003366+0.000660*k^2-0.000089*k^4;
a = Ebcc+da;
b = 0.565004-0.026134*k^2-0.002689*k^4;
c = -0.206893-0.086384*k^2+0.018278*k^4;
d = -0.031402+0.042429*k^2-0.008037^4;
u = a*gam+b*gam^(1/3)+c+d*gam^(-1/3);
UiiEQ = u/gam+k/2;
end

