function [xData,yData,yDataErr,plotInfo] = extractEnergyData(energyData,simParams,hc,flags)
%% Program notes

%   This program extracts energy data (e.g. temperatures, correlation energy, potential energy, etc) to be plotted later
%   on based on the user-specified options contained within 'flags'.

%   Inputs
%       energyData: (table) exactly as saved by 'loadEnergyInfoV2.m' - see that program for more details
%       simParams: (struct) exactly as saved by 'loadSimParamsFromFolderName.m' - see that program for more details
%       hc: (struct) hard-coded constants important for variuos calculations
%       flags: (struct) whose fields contain user-selected options

%   Outputs
%       xData: (mx1 double) column vector containing time (in either SI or dimensionless units)
%       yData: (mxn double) columns of 'yData' contain data sets to be plotted, rows within yData should correspond to
%                           rows in xData
%       yDataErr: (mxn double) same comment as yData except this array contains the standard error of the quantities
%                               found in yData
%       plotInfo: (1x1 struct) contains information for labeling the plot later on
%                       - .xlabel: (string)
%                       - .ylabel: (string)
%                       - .legend: (1xn cell) each element contains a string, should have the same number of columns as
%                                           yData
%                       - .title: (string)

%% Extract user-selected data from the 'energyData' matrix
% Preallocate space to matrices
xData = [];
yData = [];
yDataErr = [];
plotInfo = [];
plotInfo.xlabel = {};
plotInfo.legend = {};
plotInfo.ylabel = {};
plotInfo.title = {};
% Define title of plot
plotInfo.title = 'Energy vs Time';
% Time used for x-axis: choose the units 
if flags.whichTimeToUse == 1 % Use time in us
    xData = energyData.timeInUs;
    plotInfo.xlabel = 'Time (\mus)';
elseif flags.whichTimeToUse == 2 % Use dimensionless time units
    xData = energyData.timeDim;
    plotInfo.xlabel = '\omega_p_it / 2\pi';
end
% Choose data for the y-axis (could be more than one vector)
for i = 1:length(flags.whichEnergyToPlot)
    if flags.whichEnergyToPlot(i) == 1 % Plot x temp?
        yData = [yData energyData.Tx];
        yDataErr = [yDataErr energyData.TxSE];
        plotInfo.legend{length(plotInfo.legend)+1} = 'T_x';
    elseif flags.whichEnergyToPlot(i) == 2 % Plot y temp?
        yData = [yData energyData.Ty];
        yDataErr = [yDataErr energyData.TySE];
        plotInfo.legend{length(plotInfo.legend)+1} = 'T_y';
    elseif flags.whichEnergyToPlot(i) == 3 % Plot z temp?
        yData = [yData energyData.Tz];
        yDataErr = [yDataErr energyData.TzSE];
        plotInfo.legend{length(plotInfo.legend)+1} = 'T_z';
    elseif flags.whichEnergyToPlot(i) == 4
        yData = [yData energyData.Uii];
        yDataErr = [yDataErr energyData.UiiSE];
        plotInfo.legend{length(plotInfo.legend)+1} = 'U_i_i';
    elseif flags.whichEnergyToPlot(i) == 5
        yData = [yData energyData.dE];
        yDataErr = [yDataErr energyData.dESE];
        plotInfo.legend{length(plotInfo.legend)+1} = '\DeltaE';
    end
end
% Choose the units for the data on y-axis
if flags.energyUnits == 1 % Use energy units of K
    yData = convDimEnergyUnitsToTempInK(yData,hc,simParams);
    yDataErr = convDimEnergyUnitsToTempInK(yDataErr,hc,simParams);
    plotInfo.ylabel = 'T (K)';
elseif flags.energyUnits == 2 % Keep the units as dimensionless
    plotInfo.ylabel = 'Energy (E_c)';
end

end

