function [ paa ] = mcx2paa( freqs, freqResp , desiredFreq, reImFlag,WT)
%mcx2paa Returns phase and amplitude or real/imaginary parts at the requested frequencies

%   INPUTS:
% freqs         //Frequencies calculated from the MC simulations
% freqResp      //Calculated complex frequency response 
% desiredFreq   //Frequencies used in the measurement
% reImFlag      //0=phase and log(amplitude), 1 = real/imaginary parts, 2=phase
%               //and amplitude
% OUTPUTS:
% paa //Phase and amplitude (or real and imaginary parts depending on
%     //reImFlag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Interpolate to get the desired frequencies
    thisInterp = interp1(freqs, freqResp,desiredFreq,'spline',0);
    
    if length(WT) == 1
        WT = ones(length(desiredFreq)*2,1);
    end
    
    %Calculate phase and scaled amplitude
    if reImFlag == 0
        amp = abs(thisInterp);
        phase = atan2(imag(thisInterp),real(thisInterp));
        paa(1:length(desiredFreq)) = log(amp).*WT(1:length(desiredFreq)); %Scale the amplitude by natural log so it is on similar scale to phase
        paa(length(desiredFreq)+1:2*length(desiredFreq)) = phase.*WT(length(desiredFreq)+1:2*length(desiredFreq));
    %Calculate phase and unscaled amplitude
    elseif reImFlag == 2
         paa(1:length(desiredFreq))= abs(thisInterp).*WT(1:length(desiredFreq));
         paa(length(desiredFreq)+1:2*length(desiredFreq)) = atan2(imag(thisInterp),real(thisInterp)).*WT(length(desiredFreq)+1:2*length(desiredFreq));
    %Just return real and imaginary parts
    else
        paa(1:length(desiredFreq)) = real(thisInterp).*WT(1:length(desiredFreq));
        paa(length(desiredFreq)+1:2*length(desiredFreq)) = imag(thisInterp).*WT(length(desiredFreq)+1:2*length(desiredFreq));
    end
    paa = paa';
end

