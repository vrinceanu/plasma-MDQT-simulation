function [flags] = askUserForInput(directory)
%% Function Notes

    % This program will ask the user a series of questions regarding what the
    % user wants plotted. Some simulation folders can contain up to  ~100 Gb
    % of data depending on the density and the maximum simulation time,
    % therefore these questions allow the program to only have to load
    % exactly what the user wants to see. The questions are pretty
    % self-explanitory (just read the text within the 'input' functions to
    % figure out what it's asking...

    % Input:
        % - directory: (nx1 cell) whose elements are full paths to MDQT data folders, n must be greater than 0.

    % Output:
        % - flags: A structure whose fields contain user-selected program options. Some fields are optional (see below).

%% Prompt user to select options
% Initialize variable
flags.plotForSpecificV = [];
% Main option - select what type of plot to do
flags.whatToPlot = input('What do you want to plot?\n (1) Energy vs time\n (2) State Population vs Velocity for a particular time\n (3) Fit Temps to Rate Equation Model\n\n Input:');
% Depending on which plot type is selected, choose additional options.
if flags.whatToPlot == 1
    % Choose one or more energies to plot
    flags.whichEnergyToPlot = input('\nWhich energy to plot?\n (1) Tx\n (2) Ty\n (3) Tz\n (4) Uii (correlation energy)\n (5) Change in total energy\n\n Input a vector with one or more:'); 
    % Select energy units for y-axis
    flags.energyUnits = input('\nChoose energy units.\n (1) Temperature (K)\n (2) Dimensionless units of Ec=e^2/(4*pi*eps*a)\n\n Input:'); 
    % Select time units for x-axis
    flags.whichTimeToUse = input('\nChoose units for x-axis (time).\n (1) Time in us\n (2) Dimensionless units of omega_pi^-1/(2*pi)\n\n Input:');
elseif flags.whatToPlot == 2
    % Choose which time point to plot velocities for
    flags.timeForPlot = input('\nSpecify time in us for state vectors to be plotted for.\n\n Input:');
    % Choose one or more states to be plotted
    flags.whichStateToPlot = input('Choose one or more states to plot the population for.\n (1) s-state\n (2) p-state\n (3) d-state\n\n Input:');
    % Select a specific velocity to plot for if desired
    flags.plotForSpecificV = input('Plot for specific velocity (ignore by entering nothing)? Enter velocity in units m/s:\n\n Input: ');
elseif flags.whatToPlot == 3
    % Choose a window in time to be considered for fitting to rate equations
    flags.timeWindow = input('Specify time window for fitting to rate equations (e.g. [Tmin Tmax] in us)\n\n Input: ');
end
% In the event that more than one data folder is being analyzed, the legend needs to be updated to include not only the
% thing being plotted but which data set it corresponds to. 
if length(directory) > 1
    % Choose which quantity to add to the legend
    flags.baseLegendOffWhat = input('Show what quantity in the legend?\n (1) Density\n (2) Te\n (3) Initial Plasma Size\n (4) Rabi Freq for SP Trans\n\n Input:');
end

end