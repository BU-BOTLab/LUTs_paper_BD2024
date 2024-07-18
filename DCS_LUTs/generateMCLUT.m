function [LUT] = generateMCLUT(allMua, allMus, allBFI, options)
% Generate LUTs, returns structure of look up tables with values at allMua and allMus
% Also saves a copy of the LUT in the LUTs folder

%INPUTS:
% allMua        -- list of all absorption values to calculate
% allMus        -- list of all scattering values to calculate
% allBFI        -- list of all flow rates values to calculate
% rhos          -- list of all source/detector separations to calculate
% simDir        -- directory that has all the simulations
% simName       -- simualtion general name

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
saveName    = sprintf('DCSLUT_%dmm_n%d.mat',options.rho,length(allBFI));
tau         = [linspace(9.6e-7,1.7067e-06,8),...
               linspace(1.9200e-06,3.4133e-06,8),...
               linspace(3.8400e-06,6.8267e-06,8),...
               linspace(7.6800e-06,1.3653e-05,8),...
               linspace(1.5360e-05,2.5600e-05,7),...
               linspace(3.0720e-05,5.4613e-05,8),...
               linspace(6.1440e-05,1.0923e-04,8),...
               linspace(1.2288e-04,2.1845e-04,8),...
               linspace(2.4576e-04,4.3691e-04,8),...
               linspace(4.9152e-04,8.7381e-04,8)];
k0          = 2*pi*options.n/(852e-6);

%% Preallocation

% Turn vector into grid
[X,Y,Z]= meshgrid(allMua, allMus, allBFI);

% Preallocate G2 normalized
G2 = cell(size(X,1),size(X,2));

%% Generate the forward data and save it
checkMusp = 0;

for i = 1:size(X,1) % use parfor for speed
    
    
    if checkMusp == 0
        % pulls out the paths for the specific musp
        [paths, moms, lowMuspIdx, q] = pathlengthCalc(Y(i,1,1), options.muspVec, options);
        
        % calculates the G2
        G2Vec = forwardMC(X(1,:,1), Z(1,1,:), paths, moms, q, tau, k0);
        
        %updates check so this doesnt happen until new musp ranges
        checkMusp = 1;
    end
    
    % interpert G2 values if musp is between simulations
    if q == 1
        for k = 1:length(X(1,:,1))
            G2{i,k} = G2Vec{k}';
        end
        checkMusp = 0;
    else
        for k = 1:length(X(1,:,1))
            for j = 1:length(Z(1,1,:))
                G2{i,k}(j,:) = interp1([options.muspVec(lowMuspIdx),options.muspVec(lowMuspIdx+1)],[G2Vec{k,1}(:,j),G2Vec{k,2}(:,j)]',Y(i,1,1));
            end
        end
    end
    
    % checks if musp to be calculalted is out of bounds
    if i ~= size(X,1)
        if Y(i+1,1,1) >= options.muspVec(lowMuspIdx+1)
            checkMusp = 0;
        end
    end
    
end

%%  Save LUT the structure and the MAT file
LUT.G2 = G2;
LUT.rho = options.rho;
LUT.muaGrid = X;
LUT.musGrid = Y;
LUT.bfiGrid = Z;
save(fullfile(options.saveDir,saveName),'LUT','-v7.3');


