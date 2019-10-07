function [statePopulations] = loadStatePopulationsForSpecifiedTime(directory,timeStep,jobs,hc,simParams)
%% Function Notes

%% Function
numOfJobs = length(jobs);
% Convert time step to the string format of the file name
timeStepString = sprintf('%06d', timeStep);

% Define velocity bins - bin spacing should be 0.0043*a_ws*w_pi
aWS = (3/(4*pi*simParams.n))^(1/3);
wPI = sqrt(simParams.n*hc.e^2/(hc.mi*hc.eps));
binSpacing = 0.008*aWS*wPI;
velRange = [-40:binSpacing:40]/aWS/(wPI/sqrt(3));

% Loop through each job file
statePopAllJobs = cell(length(velRange),1);
for i = 1:numOfJobs % iterate through jobs
    currPopFile = dlmread([directory '\job' num2str(jobs(i)) '\statePopulationsVsVTime' timeStepString '.dat']);
    statePopForThisJob = cell(length(velRange),1);
    for j = 1:size(currPopFile,1)   % iterate through each particle for the current job
        % Find the velocity within velRange that the current velocity is
        % closest to
        currVel = currPopFile(j,1);
        velDif = abs(currVel-velRange);
        [~,minInd] = min(velDif);   % this finds the velocity bin index that best represents the current velocity
        
        % Bin velocities into currStatePop
        statePopForThisJob{minInd,1} = [statePopForThisJob{minInd,1};currPopFile(j,2:4)];
    end
    % Calculate mean population
    for k = 1:length(velRange)
        if isempty(statePopForThisJob{k,1}) == 0
            statePopAllJobs{k,1} = [statePopAllJobs{k,1}; mean(statePopForThisJob{k,1},1)];
        end
    end
    
end
% Create table of state populations
statePopulations = array2table([]);
statePopulations.vel = cell(length(velRange),1);
statePopulations.s = cell(length(velRange),1);
statePopulations.sSTD = cell(length(velRange),1);
statePopulations.sSTE = cell(length(velRange),1);
statePopulations.p = cell(length(velRange),1);
statePopulations.pSTD = cell(length(velRange),1);
statePopulations.pSTE = cell(length(velRange),1);
statePopulations.d = cell(length(velRange),1);
statePopulations.dSTD = cell(length(velRange),1);
statePopulations.dSTE = cell(length(velRange),1);

for k = 1:length(statePopAllJobs)
    statePopulations.vel{k} = velRange(k);
    if isempty(statePopAllJobs{k,1}) == 0
        currPopSize = size(statePopAllJobs{k});
        currMeanPops = mean(statePopAllJobs{k},1);
        currSTDPops = std(statePopAllJobs{k},0,1);
        currSTEPops = currSTDPops./sqrt(currPopSize(1));
        statePopulations.s{k} = currMeanPops(1);
        statePopulations.sSTD{k} = currSTDPops(1);
        statePopulations.sSTE{k} = currSTEPops(1);
        statePopulations.p{k} = currMeanPops(2);
        statePopulations.pSTD{k} = currSTDPops(2);
        statePopulations.pSTE{k} = currSTEPops(2);
        statePopulations.d{k} = currMeanPops(3);
        statePopulations.dSTD{k} = currSTDPops(3);
        statePopulations.dSTE{k} = currSTEPops(3);
    end
end

% Ensure all the rows of statePopulations actually contain populations
ind = [];
for k = 1:length(statePopulations.vel(:))
    if isempty(statePopulations.s{k})
        ind = [ind k];
    end
end

statePopulations(ind,:) = [];

end

