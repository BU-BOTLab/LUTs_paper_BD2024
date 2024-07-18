% Code to create lookup tables based on information from the MC simulations
% Written by Carlos Gomez
% Written on 7/23/21
% Last edited on 7/23/21

clear
close all
clc

%% bookkeeping
startFreq = 139;
freqStep = 10;

options.saveDir  = '..\..\Data\LUTs\DOS\Subject';
options.simDir   = '..\..\Data\Simulations\Subject';
options.fname    = 'semiInfRef';
options.muspVec  = 0.1:.1:1.5;
options.rho      = 25;
options.NA       = 0.66;
options.n        = 1.37;
options.unitInmm = 1;

%% generate LUTs

muaVec  =[.005:.0001:.01,.011:.001:.1]; 
muspVec = .1:.01:1.5;

tic
LUT730 = generateMCLUT(muaVec,muspVec,startFreq         ,options);
toc
tic
LUT850 = generateMCLUT(muaVec,muspVec,startFreq+freqStep,options);
toc
