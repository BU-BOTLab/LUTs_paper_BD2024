function [ G2Vec ] = forwardMC( muaVec, bfiGuess, paths, moms, q, tau, k0)
% Calcualtes thee tpsf based on the photon pathlenghts and photon weight.
% The TPSF is then FFT in order to get all of the freqs and the calculates
% the real and imag components of the signal at a specific freq.

%INPUTS:
% muaVec        -- vector of all the differen mua in the LUT
% bfiGuess      -- bfi in mm2/s
% paths         -- the pathlengths from the white monte carlo simulations
% moms          -- the momentum from the white monte carlos simulations
% q             -- flag to tell if there interpertation needs to occur between simulations
% tau           -- time delays for the autocorrelation
% k0            -- wavenumber

%OUTPUTS:
% G2          -- Normalized G2 curve 

%% Calculations of various photon derived variables

for i = 1:q
    
    % calculates photon weight
    photonWeight = exp(-muaVec.*paths{i});
    photonWeight2 = repmat(photonWeight',[1,1,length(tau)]);
    photonWeight2 = permute(photonWeight2,[2,3,1]);
    
    % calculate photon momentum
    photonMomentum = (moms{i} * 6 * tau) .* bfiGuess;
    photonMomentum2 = exp(-1/3*k0^2*photonMomentum);
    
    %calculate G2
    for j = 1:length(muaVec)
            G1 = squeeze(sum((photonMomentum2.*photonWeight2(:,:,j)),1));
            g1 = G1./G1(1);

            G2Vec{j,i} = g1.^2;
    end

end

end

