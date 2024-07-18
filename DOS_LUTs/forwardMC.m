function [ paa ] = forwardMC( muaVec, muspGuess, lowMuspIdx, muspVec, freqs, paths, fs, c_mmps, edges, q, reImFlag, WT)
% Calcualtes thee tpsf based on the photon pathlenghts and photon weight.
% The TPSF is then FFT in order to get all of the freqs and the calculates
% the real and imag components of the signal at a specific freq.

%INPUTS:
% muaVec        -- vector of all the differen mua in the LUT
% muaGuess      -- absorption coefficient in 1/mm
% lowMuspIdx    -- the index of the lowest simulation
% muspVex       -- vector space of musp used in simulation
% freqs         -- the desired modulation frequencies to return in MHz
% paths         -- the pathlengths from the white monte carlo simulations
% fs            -- the sampling rate
% c_mmps        -- speed of light in mm/s
% q             -- flag to tell if there interpertation needs to occur between simulations
% edges         -- the bins for the fft
% reImFlag      -- if 0 returns phase and log(amp), if 1, returns real/imag parts, if 2 returns phase and amp

%OUTPUTS:
% paa           -- Phase and amplitude of the given OPs at freqs

%% Calculations of various photon derived variables
binVal = zeros(length(edges),length(muaVec),q);


for i = 1:q

    [~,n] = size(paths{i});

    if n > 1
        paths{i}(:,1) = [];
    end

    % calculates total path length
    totalPathLength = sum(paths{i},2);
    
    % calculates photon weight
    photonWeight = exp(-muaVec.*paths{i});
    
    % calculate time the photon was detected
    times = paths{i}/c_mmps;
    
    % calcualtes total photons
    for d = 2:length(edges)
        thesePhotons = find(times > edges(d-1) & times <= edges(d));
        binVal(d-1,:,i) = sum(photonWeight(thesePhotons,:),1); % Carlos
    end
end

%Count the number of photons that arrive in each temporal bin
if q == 1
    interpTPSF = binVal;
else
    for i = 1:length(muaVec)
        interpTPSF(:,i) = interp1([muspVec(lowMuspIdx),muspVec(lowMuspIdx+1)],squeeze(binVal(:,i,:))',muspGuess);
    end
end

%Pad the interpolated TPSF and performs FFT on it
paddedTPSF = zeros(2^15,length(muaVec));
paddedTPSF(1:size(interpTPSF),:) = interpTPSF;
fResp = fft(paddedTPSF);

%Calculate frequency bins
bins = ((0:length(paddedTPSF)-1) * fs / length(paddedTPSF)) / 1e6;

%% Get the amplitude and phase at the desired frequencies
for i = 1:length(muaVec)
    [paa(:,i)] = mcx2paa(bins,fResp(:,i)',freqs,reImFlag,WT);
end

end

