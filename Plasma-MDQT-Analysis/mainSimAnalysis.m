%% Initialize Program
% Clear workspace, figures, and command window
clc
close all
clear
% Add location of subprograms to search path - make sure this is correct!
addpath('C:\Users\gg25\Documents\GitHub\Plasma-MDQT-Analysis\subPrograms')

%% Function Notes

    % This program is the user interface for plotting data saved by the MDQT simulation program. Please see the ReadMe file
    % for more information about the MDQT code and and this program.

    % All you need to do to run this program is change the user input options below and click run. The program will then 
    % ask you a series of questions about what you want to plot.

%% User Input
%%% Specify Directory to Analyze %%%
    % Note: The user may select a folder which contains one or more simulation folders (e.g. a folder that begins with 'Ge...')
    % or the user may select one or more simulation folders. In the event that more than one simulation folder is chosen, the
    % program will request that the user specify which input parameter to use in the legend to distinguish the different data
    % sets.
[directory] = uigetfile_n_dir('F:\LaserCoolingSimulationsForPaper\Grant');
useErrorBars = 1; % Plot data with (1) or without (0) standard error bars.

%% Program responsible for the plotting
%   This is the main program that loads the data and does the plotting.
[directory,f,ax] = analyzeMDQT(directory,useErrorBars);