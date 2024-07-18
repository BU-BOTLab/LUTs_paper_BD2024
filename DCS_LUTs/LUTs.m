% Code to create lookup tables based on information from the MC simulations
% Written by Carlos Gomez
% Written on 7/23/21
% Last edited on 7/23/21

clear
close all
clc

%% bookkeeping
options.saveDir  = '..\..\Data\LUTs\DCS\Subject';
options.simDir   = '..\..\Data\Simulations\Subject';
options.fname    = 'semiInfRef';
options.muspVec  = 0.1:.1:1.5;
options.rho      = 25;
options.NA       = 0.16;
options.n        = 1.37;
options.unitInmm = 1;

%% generate LUTs

muaVec  =[.005:.0001:.01,.011:.001:.1]; 
muspVec = .1:.01:1.5;
bfiVec  = [.1e-6:.001e-6:1e-6,1.01e-6:.01e-6:10e-6,10.1e-6:.5e-6:40e-6]; % SCM range

tic
LUT850 = generateMCLUT(muaVec,muspVec,bfiVec,options);
toc
