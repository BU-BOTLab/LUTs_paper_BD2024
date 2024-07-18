% Code to run MLMC for no MUA and for 2 wavelengths 730 and 850

%==========================================================================
%  Large detector  pencil beam in semi infinite medium.  This should be the
%  simplest MC case that I can compare with an analytical solution in order
%  to test diffuse reflectance calculations
%==========================================================================

close all
clc
clear
clear cfg
clear cfgVec

%% subject values

subNum = [17,15,12,18,14,19,20,13,21,22,23,24,25,26,27,28,29];
subSkin = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]; %based on ITA
subThickness = [2,3,2.3,2.9,2.2,1.5,1.4,3.8,1.8,3.1,1.8,2.7,3.2,1.5,2.2,3.3,2.2]*10;

%% Bookkeeping
addpath(genpath('/project/botlab/Carlos/mcxlab'))
saveDir = '../../Data/BU/Simulations/5618E_0';

laserToSim = [850];
musRange   = 1:1:7;    % This is mus NOT musp!

%% Build up model

% simulation paramters & settings
cfg.gpuid        = '11';
cfg.autopilot    = 1;
cfg.unitinmm     = 0.1;     % lenght of a voxel in mm
cfg.isreflect    = 1;       % Reflect index mitmatch
cfg.maxdetphoton = 1e7;     % save up to this many photons (effects memory significantly)
cfg.issaveexit   = 1;       % Save the position and direction of detected photons (can't hurt)
cfg.issaveref    = 1;       % Also save diffuse reflectance
cfg.ismomentum   = 1;       % Also save mometum transfers
cfg.nphoton      = 1e9;     % number of photons

% source and detector
cfg.srctype   = 'gaussian';
cfg.srcpos    = [180, 300, 1];           % source position
cfg.srcdir    = [000, 000, 1];           % source direction
cfg.srcparam1 = [35/2, 0, 0, 0];
cfg.detpos    = [430, 300, 2, 35/2];     % detector position and radius (x,y,z,r)

% Timing information
cfg.tstart = 0;       % Time the simulation starts
cfg.tstep  = 5e-9;    % Steps to take
cfg.tend   = 5e-9;    % When to end.  The output will have [tstart:tstep:tend] slices at each of the different time points

for w = 11
    

    %% model creation
    cfg.vol             = uint8(ones(600,600,1000)) * 3;     % creates volume and sets it all to fourth layer (muscle)
    cfg.vol(:,:,1)      = 0;                                 % pad a layer of 0s to get diffuse reflectance
    cfg.vol(:,:,2:11)   = 1;                                 % skin: 1 mm layer
    cfg.vol(:,:,12:11+subThickness(w))  = 2;                 % lipid layer


    cfg.prop = [ 0,  0,  1,  1;      % air
                 0,  0,  0,  0;      % Epidermis/dermis will vary per wavelength
                 0,  0,  0,  0;      % Lipid will vary per wavelength
                 0,  0,  0,  0];     % Muscle (will vary)[mua,mus,g,n]

    %% Model: creates a stucutre with multiple models to run (different wavelengths and updated mua & mus)

    for SkinTone = subSkin(w)
        iter = 0;
        for k = laserToSim
            for j = 1:length(musRange)
                iter = iter + 1;

                cfgVec(iter) = cfg;

                switch k
                    case 730

                        switch SkinTone
                            case 1
                                cfgVec(iter).prop(2,:) = [0.0753, 1.3352/.1, 0.9, 1.4];    % Epidermis/dermis @ 730 (Lowest Skintone)
                            case 2
                                cfgVec(iter).prop(2,:) = [0.0510, 1.3352/.1, 0.9, 1.4];    % Epidermis/dermis @ 730 (Medium Skintone)
                            case 3
                                cfgVec(iter).prop(2,:) = [0.1494, 1.3352/.1, 0.9, 1.4];    % Epidermis/dermis @ 730 (Highest Skintone)
                        end

                        cfgVec(iter).prop(3,:) = [0.01, 1.19/.1, 0.9, 1.44];            % Lipid @ 730

                    case 830

                        switch SkinTone
                            case 1
                                cfgVec(iter).prop(2,:) = [0.0231, 1.1706/.1, 0.9, 1.4];    % Epidermis/dermis @ 830 (Lowest Skintone)
                            case 2
                                cfgVec(iter).prop(2,:) = [0.0172, 1.1706/.1, 0.9, 1.4];    % Epidermis/dermis @ 830 (Medium Skintone)
                            case 3
                                cfgVec(iter).prop(2,:) = [0.0435, 1.1706/.1, 0.9, 1.4];    % Epidermis/dermis @ 830 (Highest Skintone)
                        end

                        cfgVec(iter).prop(3,:) = [0.01, 1.09/.1, 0.9, 1.44];             % Lipid @ 830
                    case 850

                        switch SkinTone
                            case 1
                                cfgVec(iter).prop(2,:) = [0.0126, 1.1424/.1, 0.9, 1.4];    % Epidermis/dermis @ 850 (Lowest Skintone)
                            case 2
                                cfgVec(iter).prop(2,:) = [0.0105, 1.1424/.1, 0.9, 1.4];    % Epidermis/dermis @ 850 (Medium Skintone)
                            case 3
                                cfgVec(iter).prop(2,:) = [0.0224, 1.1424/.1, 0.9, 1.4];    % Epidermis/dermis @ 850 (Highest Skintone)
                        end

                        cfgVec(iter).prop(3,:) = [0.01, 1.07/.1, 0.9, 1.44];             % Lipid @ 850

                end

                cfgVec(iter).prop(4,:) = [0, musRange(j), 0.9, 1.37]; %This is to replicate the work by Fang 2009
            end
        end


        %% Model: simulation and save data
        iter = 0;

        for i = 1:length(laserToSim)
            for k = 1:length(musRange)
                iter = iter + 1;

                %simulation command
                [flux,det] = mcxlab(cfgVec(iter));

                %saving data
                thisFname  = sprintf('semiInfRef_%03dnm_%02d_of_%02d.mat',laserToSim(i),k,length(musRange));
                save(fullfile([saveDir,num2str(subNum(w))],thisFname),'det','cfgVec','-v7.3');

            end
        end
    end

end

%% more bookkeeping
rmpath(genpath('B:\Users\Carlos\Projects\MonteCarlo\mcxlab'))
