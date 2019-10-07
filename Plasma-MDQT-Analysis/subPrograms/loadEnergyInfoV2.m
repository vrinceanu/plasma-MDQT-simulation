function [data] = loadEnergyInfo2(directory,jobs)
%% FUNCTION NOTES
    % This function extracts all the information from the 'energies.dat' files within each job folder. Each 
    % 'energies.dat' file contains 7 columns of information where each row represents a different time within the
    % simulation. The columns are organized as: [t Ex Ey Ez Uii dEtot], where t is time, E is kinetic energy along the
    % specified axis, Uii is total potential energy, and Etot is the total energy. All energies represent the energy per
    % particle (e.g. we add up the x-kinetic energy for all particles and divide by total number of particles).
    
    % Units within 'energies.dat' file
        % time is recorded in units of wpE^-1 = sqrt(n*e^2/(3*m*epsilon))
        % energy is recorded in units of Ec = e^2/(4*pi*eps*a), where a is the wigner-seitz radius
        
    % Outputs
        % This function contains a single output - 'data' - which is a table containing the loaded energy vs time information
        % that is averaged over all jobs. We also provide the standard deviation and standard error for each quantity (other
        % than time) within the table. The structure of the table should be self-explanitory after reviewing the 'Create Data
        % Table with Statistics' section below - in particular pay attention to the units.

%% LOAD ENERGY DATA
%   Read in the 'energies.dat' file in each job folder
%       - no calculations or conversions are done in this section
%       - this section simply reads in data from the matrices
%       - final output of this step is 'energyInfo'
%       - see next section for description of what each column entails
numOfJobs = length(jobs);

%%% Import data from each 'energies.dat' file for each job
loadedData{length(jobs)} = [];
for i = 1:length(jobs)
    currEnergyMatrix = dlmread([directory '\job' num2str(jobs(i)) '\energies.dat']);
    loadedData{i} = currEnergyMatrix;
end

%%% Ensure that number of timesteps (number of rows) for each
%%% 'energies.dat' matrix is the same. For some reason they differ by one sometimes...not sure why
clear numTimeSteps
numTimeSteps = zeros(size(jobs));
for i = 1:length(jobs)
    numTimeSteps(i) = size(loadedData{i},1);
end
numTimeSteps = min(numTimeSteps); % Identify minimum number of time steps

%   Keep only the rows of 'energies.dat' up to numTimeSteps 
energyInfo = [];
energyInfo(numTimeSteps,7,length(jobs)) = 0;
for i = 1:length(jobs)
    currEnergyMatrix = loadedData{i};
    energyInfo(1:numTimeSteps,:,i) = currEnergyMatrix(1:numTimeSteps,:);
end

%% CREATE DATA TABLE WITH STATISTICS
% For each column of the 'energies.dat' file, calculate the mean, std, and
% ste.
data = array2table([]);

% This is time in units of w_pE^-1, where w_pE = w_pi/sqrt(3) = sqrt(n*e^2/(3*m*eps))
data.time = mean(energyInfo(:,1,:),3);

% Import kinetic energy - note what I actually record is the temperature
    % For kinetic energy in units of Ec and temperature in units of Ec/kB, the relation between the two quantities is
    % KE = 0.5*T. The 'energies.dat' matrix records the kinetic energy. I am interested in recording the temperature, so I
    % multiply columns 2-4 of the 'energies.dat' matrix by 2.
    
    % The temperatures recorded below are in units of Ec/kB
    
data.Tx = mean(energyInfo(:,2,:),3).*2;
data.TxSD = std(energyInfo(:,2,:),0,3).*2;
data.TxSE = data.TxSD./sqrt(numOfJobs);

data.Ty = mean(energyInfo(:,3,:),3).*2;
data.TySD = std(energyInfo(:,3,:),0,3).*2;
data.TySE = data.TySD./sqrt(numOfJobs);

data.Tz = mean(energyInfo(:,4,:),3).*2;
data.TzSD = std(energyInfo(:,4,:),0,3).*2;
data.TzSE = data.TzSD./sqrt(numOfJobs);

% Record correlation energy normalized by particle number in units Ec
    % This was discussed in the MDQT paper, but if U is the total potential energy per particle, then Uii(t) = U(t)-U(0)
    % represents the correlation energy per particle (assuming that U(0) represents a completely uncorrelated gas - which in
    % our case it does because particle positions are initially randomized).
Ki = energyInfo(1,2,:)+energyInfo(1,3,:)+energyInfo(1,4,:); % this is the initial kinetic energy, which is zero for normal simulations
Uii = energyInfo(:,5,:) - (energyInfo(1,5,:)+Ki);
data.Uii = mean(Uii(:,:,:),3);
data.UiiSD = std(Uii(:,:,:),0,3);
data.UiiSE = data.UiiSD./sqrt(numOfJobs);

% This is the difference between the total energy (Potential+Thermal) at time t and
% the initial total energy. This quantity becomes increasingly negative as
% t increases during the laser cooling process because energy is being
% removed. This quantity is in units of Ec = e^2/(4*Pi*Eps*a)
K = energyInfo(:,4,:)+energyInfo(:,2,:)+energyInfo(:,3,:);
dE = K+energyInfo(:,5,:)-Ki-energyInfo(1,5,:);
data.dE = mean(dE,3);
data.dESD = std(dE,0,3);
data.dESE = data.dESD./sqrt(numOfJobs);

% data.dE = mean(energyInfo(:,6,:),3);
% data.dESD = std(energyInfo(:,6,:),0,3);
% data.dESE = data.dESD./sqrt(numOfJobs);

% This is the expansion velocity of the plasma. This can be ignored if
% you're considering the center of the plasma, but is relevant if you're
% considering an expanding plasma. See the 'loadSimParamsFromFolderName'
% function for more details on that. This quantity is in units of
% a*w_pE
data.v = mean(energyInfo(:,7,:),3);
data.vSD = std(energyInfo(:,7,:),0,3);
data.vSE = data.vSD./sqrt(numOfJobs);

end

