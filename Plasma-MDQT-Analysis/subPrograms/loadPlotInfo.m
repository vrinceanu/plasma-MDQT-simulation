function [xData,yData,yDataErr,plotInfo] = loadPlotInfo(directory,jobs,simParams,hc,flags)
%% Function Notes

    % This function is responsible for loading whatever the user chose to plot. The plotting is not actually done here.
    
    % Inputs
        % directory: (string) full path to a single simulation data folder
        % jobs: (nx1 double) each element contains a unique job number, where n is the total number of jobs
        % simParams: (struct) each field contains a simulation input parameter in SI units
        % hc: (struct) each field contains a hard-coded symbolic constant for use in various calculations in SI units
        % flags: (struct) fields contain user-selected program options
        
    % Outputs
        % xData: (nx1 double) contains data for x-axis of plot
        % yData: (nxm double) contains corresponding data m data sets to be plotted on y-axis
        % yDataErr: (nxm double) contains standard error of data to be plotted on y-axis
        % plotInfo: (struct) fields contain plot label information

%% Load energy data
% The time vector is only contained within the 'energies.dat' file, so I'm
% going to load that first. 
[energyData] = loadEnergyInfoV2(directory,jobs);

%% Convert time units
% Calculate real time - the simulation time (as recorded within the data
% table) is in units of omega_pE^-1. This converts time to units of us.
energyData.timeInUs = energyData.time(:).*sqrt(3*hc.mi*hc.eps/(simParams.n*hc.e^2))*1e6;

% Calculate dimensionless time for plotting (we want to plot
% omega_pi*t/2pi). This converts the program time to that
energyData.timeDim = sqrt(3)/(2*pi).*energyData.time(:);

%% Load requested data into xData, yData, and plotInfo
% Initialize matrices
xData = [];
yData = [];
yDataErr = [];
plotInfo = [];
% Load data into matrices
if flags.whatToPlot == 1
    [xData,yData,yDataErr,plotInfo] = extractEnergyData(energyData,simParams,hc,flags);
elseif flags.whatToPlot == 2
    [xData,yData,yDataErr,plotInfo] = extractStatePopForParticularT(directory,jobs,energyData.timeInUs,hc,simParams,flags);
elseif flags.whatToPlot == 3
    [xData,yData,yDataErr,plotInfo] = fitTempsToRateEqns(energyData,hc,simParams,flags);
end

end

