function [xData,yData,yDataErr,plotInfo] = extractStatePopForParticularT(directory,jobs,timeInUs,hc,simParams,flags)
%% Program notes

%   This program loads the state population data as requested by the user.

%   Inputs
%       directory: (string) full path to a MDQT data folder
%       jobs: (nx1 double) vector whose elements contain each unique job number, n is total number of jobs
%       timeInUs: (1x1 double) user-specified time for which to load state populations for
%       hc: (1x1 struct) hard-coded quantities useful for some calculations - see 'analyzeMDQT.m'
%       simParams: (1x1 struct) initial simulation parameters loaded from 'directory' by 'loadSimParamsFromFolderName.m'
%       flags: (1x1 struct) user-selected program options

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

%% Load state populations from 'directory'
%%% Preallocate space to matrices
yData = [];
yDataErr = [];
plotInfo = [];
plotInfo.xlabel = {};
plotInfo.legend = {};
plotInfo.ylabel = {};
plotInfo.title = {};

%%% Obtain time step closest to specified time
[timeStep,timeForPlotInUs] = timeStepClosestToSpecifiedTime(flags.timeForPlot,timeInUs);

%%% Load the State Populations for the specified time
[statePopulations] = loadStatePopulationsForSpecifiedTime(directory,timeStep,jobs,hc,simParams);

%%% Use velocity along x-axis, convert to SI units
xData = convDimVelToSI(cell2mat(statePopulations.vel(:)),hc,simParams);
plotInfo.xlabel = 'Velocity (m/s)';

%%% Specify label for y-axis
plotInfo.ylabel = 'State Population';

%%% Specify plot title
plotInfo.title = ['State Pop vs Vel for t = ' sprintf('%.2f',timeForPlotInUs) ' \mus'];

%%% Choose which state populations to plot
for i = 1:length(flags.whichStateToPlot)
    if flags.whichStateToPlot(i) == 1
        yData = [yData cell2mat(statePopulations.s)];
        yDataErr = [yDataErr cell2mat(statePopulations.sSTE)];
        plotInfo.legend{length(plotInfo.legend)+1} = 'P_s(v)';
    elseif flags.whichStateToPlot(i) == 2
        yData = [yData cell2mat(statePopulations.p)];
        yDataErr = [yDataErr cell2mat(statePopulations.pSTE)];
        plotInfo.legend{length(plotInfo.legend)+1} = 'P_p(v)';
    elseif flags.whichStateToPlot(i) == 3
        yData = [yData cell2mat(statePopulations.d)];
        yDataErr = [yDataErr cell2mat(statePopulations.dSTE)];
        plotInfo.legend{length(plotInfo.legend)+1} = 'P_d(v)';
    end
end

end

