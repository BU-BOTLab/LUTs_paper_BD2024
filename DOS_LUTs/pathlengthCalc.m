function [paths, lowMuspIdx, q] = pathlengthCalc(muspGuess, muspVec, options)
%% Calcualtes the photon paths from the simulation based on muspGuess

%INPUTS:
% muspGuess     -- reduced scattering coefficient in 1/mm
% muspVec       -- range of simulated reduced scattering coefficients in 1/mm
% options       -- various user defined inputs

%OUTPUTS:
% paths         -- the photon paths from the simulation based on MCX simulation

%% error checking
if muspGuess > max(muspVec)
    error('Scattering guess is higher than simulation range')
end
if muspGuess < min(muspVec)
    error('Scattering guess is lower than simulation range')
end

%% Loads simulation

lowMuspIdx = find(muspVec <= muspGuess,1,'last');

if muspGuess == muspVec(lowMuspIdx)
    numSimsToLoad = 1;
else
    numSimsToLoad = 2;
end

%% Calculates the paths
for q = 1:numSimsToLoad
    
    %loads simulation
    thisFile =sprintf('%s_%02d_of_%02d.mat',options.fname,lowMuspIdx+(q-1),length(muspVec));
    load(fullfile(options.simDir,thisFile),'det');
    
    th = asin(options.NA/options.n(1));
    z = acos(-det.v(:,3));
    paths{q} = det.ppath(z<= 2*th,:)*options.unitInmm;
    moms{q} = det.mom(z<= 2*th,:);
end

end