function [xData,yData,yDataErr,plotInfo] = fitTempsToRateEqns(energyData,hc,simParams,flags)
%% Program notes

%   This program fits the temperatures parallel/perpendicular to the laser cooling axis and the correlation energy to a
%   rate equation model described within the MDQT publication.

%   Inputs
%       energyData: (table) as loaded by 'loadEnergyInfoV2.m'
%       hc: (1x1 struct) hard-coded quantities useful for calculations - see 'analyzeMDQT.m'
%       simParams: (1x1 struct) initial simulation parameters loaded by 'loadSimParamsFromFolderName.m'
%       flags: (1x1 struct) user-selected program options - see 'askUserForInput.m'

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

%% Obtain time vector from 'energyData' within user-specified time window for fitting rate equations
% Extract time vector from 'energyData' - use time vector that is in units us
timeInUs = energyData.timeInUs;
% Prompt user to select time window - user specifies time in us
timeWindow = flags.timeWindow;
% Check to make sure user properly specified 'timeWindow'
if length(timeWindow)~=2 || timeWindow(1)>timeWindow(2) || timeWindow(1)<0 
    error('Time window not specified properly. Should be 1x2 double with times in us.');
end
% Find time within 'timeInUs' that is closest to lower time window
[minInd] = timeStepClosestToSpecifiedTime(timeWindow(1),timeInUs);
minInd = minInd+1; % I do this because the index that comes from the previous function is zero based
% Find time within 'timeInUs' that is closest to upper time window
[maxInd] = timeStepClosestToSpecifiedTime(timeWindow(2),timeInUs);
maxInd = maxInd+1; % I do this because the index that comes from the previous function is zero based
% Only keep times within the specified window
timeForFit = timeInUs(minInd:maxInd);

%% Prepare time and temperature matrices for fitting with 'nlinfit'
% Create temperature vector for fit [Tx Ty Uii/k]
tempForFit = [energyData.Tx (energyData.Ty+energyData.Tz)/2 energyData.Uii];
tempErrForFit = [energyData.TxSE (energyData.TySE+energyData.TzSE)/2 energyData.UiiSE];
% Only keep values within specified time window
tempForFit = tempForFit(minInd:maxInd,:);
tempErrForFit = tempErrForFit(minInd:maxInd,:);

%% Create initial guesses for fit
% Calculate plasma frequency in SI units
wpi = sqrt(simParams.n*hc.e^2/hc.eps/hc.mi);
% Initial guess for laser cooling rate
Beta = .03; % units us^-1
% Initial guess for cross axis thermalization
Nu = .1*wpi*1e-6; % units us^-1
% Initial guess for correlation temperature relaxation
Mu = .1*wpi*1e-6; % units us^-1
% Put it all together
beta0 = [Beta Nu Mu];

%% Create weights for fit from standard error
% Here, I'm going to use the standard error in the temperature to weight the observations in the non-linear regression
% below. See the 'nlinfit' help function for more informatioin on how weights affect the fit coefficients. From that
% information, it appears that the weights should sum to N (the number of observations) and the higher the weight the
% more influence it has over the fit. Thus, I will use the inverse of the standard error to weight each measurement.
w = 1./tempErrForFit;
w = w./sum(w).*length(w);

%% Fit the temperatures to the rate equations
% Specify initial temperature
T0 = tempForFit(1,:);
% Initialize fit model
modelfun = @(c,times) solveRateEquations(c,times,simParams.k,T0);
% Fit for the best coeffs
[coeff,R,~,CovB,~,~] = nlinfit(timeForFit,tempForFit(:),modelfun,beta0,'Weights',w(:));
ci = nlparci(coeff,R,'covar',CovB);
[~,delta] = nlpredci(modelfun,timeForFit,coeff,R,'Covar',CovB);
% Extract temperatures from fit
tempFromFit = reshape(modelfun(coeff,timeForFit),[],3);
tempErrFromFit = reshape(delta,[],3);

%% Output fitted temperatures
xData = repmat(timeForFit,1,6);
yData = [tempForFit(:,1) tempFromFit(:,1) tempForFit(:,2) tempFromFit(:,2) -tempForFit(:,3) -tempFromFit(:,3)];
yData = convDimEnergyUnitsToTempInK(yData,hc,simParams);
yDataErr = [tempErrForFit(:,1) tempErrFromFit(:,1) tempErrForFit(:,2) tempErrFromFit(:,2) -tempErrForFit(:,3) -tempErrFromFit(:,3)];
yDataErr = convDimEnergyUnitsToTempInK(yDataErr,hc,simParams);

%% Specify plot information
% Specify xLabel
plotInfo.xlabel = 'time (\mus)';
% Specify yLabel
plotInfo.ylabel = 'T (K)';
% Specify title
plotInfo.title = ['Fit Params: \beta = ' num2str(round(coeff(1),3)) ' \mus^-^1, \nu = ' num2str(round(coeff(2)/wpi*1e6,3)) '\omega_p_i, \mu = ' num2str(round(coeff(3)/wpi*1e6,3)) '\omega_p_i'];
% Specify legend
plotInfo.legend = {'T_|_|', 'T_|_| Fit', 'T_\perp', 'T_\perp Fit', '-U_i_i/k_B','-U_i_i/k_B Fit'};
% Export coefficients and their error
coeffErr = (ci(:,2)-ci(:,1))./2;
plotInfo.coeff = [coeff(1) coeff(2)/wpi*1e6 coeff(3)/wpi*1e6];
plotInfo.coeffErr = [coeffErr(1) coeffErr(2)/wpi*1e6 coeffErr(3)/wpi*1e6];

end

