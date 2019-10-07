clc
clear all
close all

filename = 'F:\LaserCoolingSimulationsForPaper\Grant\03.21.19 Investigate Use of Timestep and Sample Freq\simDarkStatePop_t7us_timestep001\Ge10Density22E+11Sig040Te19SigFrac0DetSP-100DetDP0OmSP100OmDP100NumIons3500';
numOfJobs = 99;

for i = 1:numOfJobs
    currFileName = [filename '\job' num2str(i) '\energies.dat'];
    currEnergyMat = dlmread(currFileName);
    cutEnergyMat = currEnergyMat(1:83,:);
    dlmwrite(currFileName,cutEnergyMat,'\t')
end

