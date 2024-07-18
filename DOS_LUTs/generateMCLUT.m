function [LUT] = generateMCLUT(allMua, allMus, freqs, options)
% Generate LUTs, returns structure of look up tables with values at allMua and allMus
% Also saves a copy of the LUT in the LUTs folder 

%INPUTS:
% allMua        -- list of all absorption values to calculate
% allMus        -- list of all scattering values to calculate
% freqs         -- frequency to calculate
% options       -- various user defined inputs

%OUTPUTS:
% LUT           -- Structure of LUT with following fields
% rp            -- real part of the data (size length(allMua)xlength(allMus)xnumFreqs)
% ip            -- imaginary part of the data (same size as above)
% rho           -- source-detector separation
% freqs         -- modulation frequencies
% muaGrid       -- Grid of mua values in the LUT
% musGrid       -- Grid of mus values in the LUT

disp('Generating new LUT')

%% hard coded variables 
saveName    = sprintf('DOSLUT_%dmm_%dMHz_n%d.mat',options.rho,freqs,length(allMua));
reImFlag    = 1;
endTime     = 5e-9;
secPerSamp  = 1e-11;
fs          = 1/secPerSamp;
c_mmps      = 3e11./options.n; 
edges       = 0:secPerSamp:endTime; %How many bins to use

%% Preallocation 

% Turn vector into grid
[X,Y]= meshgrid(allMua, allMus);

% Separate the real and imaginary parts
rp = zeros(size(X,1),size(X,2));
ip = zeros(size(X,1),size(X,2));

%% Generate the forward data and save it

parfor i = 1:size(X,1) % use parfor for speed
    
    % pulls out the paths for the specific musp 
    [paths, lowMuspIdx, q] = pathlengthCalc(Y(i,1), options.muspVec, options);
    
    % calculates the real and imaginary data  
    rawDat   = forwardMC(X(i,:), Y(i,1), lowMuspIdx, options.muspVec, freqs, paths, fs, c_mmps, edges, q, reImFlag,0);
    rp(i,:)  = rawDat(1,:);
    ip(i,:)  = rawDat(2,:);
end

%%  Save LUT the structure and the MAT file
LUT.rp = rp;
LUT.ip = ip;
LUT.freqs = freqs;
LUT.rho = options.rho;
LUT.muaGrid = X;
LUT.musGrid = Y;
save(fullfile(options.saveDir,saveName),'LUT','-v7.3');


