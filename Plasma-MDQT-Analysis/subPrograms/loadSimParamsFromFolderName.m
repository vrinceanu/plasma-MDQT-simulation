function [simParams] = loadSimParamsFromFolderName(directory,hc)
%% FUNCTION NOTES

    % This function extracts the simulation input parameters for the laser cooling
    % simulation code that were used to produce the data within 'directory' (should be a string).
    % The information for all these parameters is contained within the folder
    % name of 'directory'. See the ReadMe for more information.

    % The input parameter 'hc' is a structure containing hard-coded quantities
    % from the main program in SI units.

    % The output parameter 'simParams' is in SI units - all frequencies have
    % been converted to Hz. Note that \gamma is in units s^-1, technically.

%% Electron Coulomb Coupling Parameter
% The value in the folder name indicates the electron Gamma*1000
simParams.Ge = str2double(extractAfter(extractBefore(directory,'Density'),'\Ge'))/1000;

%% Ion/Electron Density
% The value in the folder name indicates the plasma density in 1x10^11 m^-3
simParams.n = str2double(extractBefore(extractAfter(directory,'Density'),'E+'))*1e11;

%% RMS Radius of the Plasma - aka the plasma size
% The value in the folder name indicates sigma*10 in mm
simParams.sigma = str2double(extractBefore(extractAfter(directory,'Sig'),'Te'))/10*1e-3;

%% Electron Temperature
% The value in the folder name is 10*Te in units K
simParams.Te = str2double(extractBefore(extractAfter(directory,'Te'),'SigFrac'))/10;

%% Location of plasma that simulation data was saved for
% The simulation code has the option to either record simulation results at
% the center of the plasma (in which case this parameter is 0) or it can
% record the laser cooling results at the location of the rms plasma size.
% If this parameter is non-zero, the program simulated laser cooling for an
% expanding cloud by calculating what sigma(t) is at every time step and
% adds an expansion velocity that obeys v_exp(r,t) = t*r/(t^2+tau^2).

% This value indicates the the fraction of the plasma size for which the
% laser cooling results were recorded for multiplied  by 100. I've saved
% this parameter as the actual fraction not multiplied by 100.
simParams.fracOfSigma = str2double(extractBefore(extractAfter(directory,'SigFrac'),'DetSP'))/100;

%% Detuning of the 408 nm Laser "Laser Cooling" Laser
% This is the detuning of the laser that addresses the S->P transition (aka
% the laser cooling laser) normalized by the natural linewidth of the
% 2S_(1/2)-2P_(3/2) transition. 

% The value in the folder name indicates the normalized DetSP*100. I the
% real value in units of Hz.
simParams.detSP = str2double(extractBefore(extractAfter(directory,'DetSP'),'DetDP'))/100*hc.gammaSP;

%% Detuning of 1033 nm and 1092 nm "Repump" Lasers
% This is the detuning of the repump lasers normalized by the natural
% linewidth of the 2S_(1/2)->2P_(3/2) transition. Note: this is not a
% mistake! This quantity is not normalized by the decay rate to the D
% state. The normalization is set by the SP transition because that is the
% relevant simulation timescale that normalizes these quantities for a
% particular simulation step.

% The value in the folder name indicates the normalized DetDP*100. I record
% the real value in Hz.
simParams.detDP = str2double(extractBefore(extractAfter(directory,'DetDP'),'OmSP'))/100*hc.gammaSP;

%% Rabi Frequency for 408 nm Laser
% This is the Rabi frequency for the 408 nm laser. As with the previous two
% quantities, this is in units of gamma_SP. The value in the folder is
% omega/gamma*100. I record the real rabi frequency in Hz.
simParams.omegaSP = str2double(extractBefore(extractAfter(directory,'OmSP'),'OmDP'))/100*hc.gammaSP;

%% Rabi Frequency of 1033 nm and 1092 nm Laser
% Same comment as the above section.
simParams.omegaDP = str2double(extractBefore(extractAfter(directory,'OmDP'),'NumIons'))/100*hc.gammaSP;

%% Number of Ions
% This is simply the number of ions the simulation used.
simParams.numIons = str2double(extractAfter(directory,'NumIons'));

%% Calculate plasma screening parameter
% This parameter is calculated based on the 'Ge' parameter contained within the folder name.
simParams.k = sqrt(simParams.Ge*3);
end

