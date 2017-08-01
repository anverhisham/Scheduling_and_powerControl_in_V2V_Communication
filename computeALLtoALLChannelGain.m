%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot for Wrapper function for link-removal algorithm 
%
% -INPUTS
%   1. 'locationPerUE' is a real vector, indicating location on a convoy. (2D location not supported!)
%   2. 'config.configChannel.isNullifySelfChannel' = false  => self-interference channel would be too high, that full duplex would be impossible 
%
% -NOTES:
%    1. Pathloss exponent n=1.9 (D5), pathLossScalingFactor=0.0025, shadowVariance = 1 (section 5.1.2 in D13) from Erik Steinmeitz Licentiate thesis
%           Erik got these values from paper 'Path Loss Modeling for Vehicle to Vehicle Communications' by Johan Karedal, Nicolai Czink
%    2. LTE Channel models (36.814: Table A.2.1.1.2-8, page:75): pathLossScalingFactor=10^-1.62, n=3.76
%
% -REFERENCES:
%   1. [1] Karedal "Pathloss modeling for vehicle-to-vehicle communications",  IEEE Trans. Veh. Technol.,  2011
%
% -\Author: Anver <anver@chalmers.se>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Assumption: 'locationPerUE' is a vector containing only real numbers!
function channelGainPerRxPerTx = computeALLtoALLChannelGain(locationPerUE,config)

    %% -Channel Parameters
    if(gf.isempty('config')) config = struct(); end
    config = patchConfig(config);
    N = numel(locationPerUE);
    
    %% -Generate Channel Gains b/w any 2 UEs, with penetration loss
    rxUEIDPerRxUEPerTxUE = [1:N]'*ones(1,N);
    txUEIDPerRxUEPerTxUE = ones(N,1)*[1:N];
    distancePerRxUEPerTxUE = nan(size(rxUEIDPerRxUEPerTxUE));
    if(~config.isWrapping)
        distancePerRxUEPerTxUE = abs(locationPerUE(rxUEIDPerRxUEPerTxUE) - locationPerUE(txUEIDPerRxUEPerTxUE));
    else    %% -Compute 'distancePerRxUEPerTxUE' after wrapping!  
        % Warning: Wrapping looses the symmetry of channel-matrix, so don't do it for Power-Control algorithms (Use only for measurements in equal-power to avoid edge effect!)
        % NOTE: Idea is to split the convoy and attach to other end, so that UE of our interest will be at index 'desiredIndexOfUE'. (This approach is adopted sice cyclic rotation will lead loss of one segment)
        desiredIndexOfUE = ceil(N/2);  %% -Every UE assumes that it's position in the convoy is this value. NOTE: Index always start from '1' not '0'
        if(any(diff(locationPerUE)<0)) gf.saveAllVariables(); gf.Error('''locationPerUE'' must non decreasing!!'); end
        totalRange = max(locationPerUE)+(locationPerUE(2)-locationPerUE(1));     %%+eps(2*max(locationPerUE));
        for iRxUE = 1:N
            locationPerUEnew = locationPerUE;
            nRightMovement = (desiredIndexOfUE-iRxUE);
            if(nRightMovement>0)        %% -Take the far right segment, flip it and attach to the left of first-UE
                locationPerUEnew(end-nRightMovement+1:end) = locationPerUE(1) - (locationPerUE(end-nRightMovement+1:end) - locationPerUE(end-nRightMovement));
            elseif(nRightMovement<0)    %% -Take the far left segment, flip it and attach to the right of last-UE  (last-UE => far right UE)
                nLeftMovement = -nRightMovement;
                locationPerUEnew(1:nLeftMovement) = locationPerUE(end) + (locationPerUE(nLeftMovement+1) - locationPerUE(1:nLeftMovement));                
            end
% %             UEatZeroLocation = find(mod([0:N-1]+(desiredIndexOfUE-iRxUE),N)==0);                %% -Find the index of UE at zero location after cyclical rotation
% %             locationPerUEnew = mod(locationPerUE-locationPerUE(UEatZeroLocation),totalRange);   %% -New locations per UE
            distancePerRxUEPerTxUE(iRxUE,:) = abs(locationPerUEnew-locationPerUEnew(iRxUE))';
        end
        distancePerRxUEPerTxUE = distancePerRxUEPerTxUE + diag(Inf*ones(1,N));
    end
    distancePerRxUEPerTxUE(distancePerRxUEPerTxUE<config.minDistanceForPathlossComputation) = config.minDistanceForPathlossComputation;     %% Lowerbounding to minimum distance for the purpose of computing pathloss
    channelGainPerRxPerTx = config.pathLossScalingFactor*distancePerRxUEPerTxUE.^-config.pathlossExponent;                      %% - Pathloss Computatin
    penetrationLossPerRxPerTx = db2pow(fix(abs(rxUEIDPerRxUEPerTxUE-txUEIDPerRxUEPerTxUE)-0.5)*config.penetrationLossInDB);
    channelGainPerRxPerTx = channelGainPerRxPerTx./penetrationLossPerRxPerTx;

    %% -Generate Shadow matrix
    if(config.isShadow)
       if(config.isRandSeed), rng('shuffle'); end
%        shadowMatrix = 10.^(shadowVariance*randn(size(channelGainPerRxPerTx)));
       shadowMatrix = lognrnd(0,config.shadowStandardDeviation,size(channelGainPerRxPerTx));
       shadowMatrix = tril(shadowMatrix,-1) + tril(shadowMatrix)';  %% -Copying Lower-Triangular part to Upper-Triangular part, to symmetricise
       channelGainPerRxPerTx = shadowMatrix.*channelGainPerRxPerTx;
    end
    if(config.isNullifySelfChannel)
        channelGainPerRxPerTx(logical(eye(size(channelGainPerRxPerTx)))) = 0;    %% Make zero self interference
    else
        channelGainPerRxPerTx(logical(eye(size(channelGainPerRxPerTx)))) = 1;    %% Make self interferring channel high
    end    
end



function config = patchConfig(config)
    configDefault.isNullifySelfChannel = true;
    %% -Karedal model Highway Scenario [1]
    configDefault.pathLossScalingFactor = 2.7542e-05;
    configDefault.pathlossExponent = 1.77;
    configDefault.shadowStandardDeviation = 3.1;        %% log(2016/12/14): Changed '3.1^2' to '3.1'. (why was the bug previously??)
    configDefault.isShadow = false;
    configDefault.isRandSeed = false;                   %% If true, then set random seed for shadow matrix
    configDefault.isWrapping = false;
    configDefault.penetrationLossInDB = 10;             %% -10dB Penetration Loss per blocking vehicle,,
    configDefault.minDistanceForPathlossComputation = 10;             %% 10m for Karedal model, 3m for Winner+B1 model as used by 3GPP 36.885 A.1.4
    %% -Merging Inputted config with Default inputs,
    config = gf.mergeTwoStructures(configDefault,config);
end


