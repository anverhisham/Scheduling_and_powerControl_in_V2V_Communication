%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% -DETAILS
%   1. Look for '#Algorithm' for knowing about algorithm.
%
% -PREREQUISITES
%   1. config must contain nonempty fields {N,gammaT,noiseVariance,Pmax}
%
% -NOTES
%   1. config.kangWangConfig.Ct = 1   => KangWang algorithm. (Otherwise it's heuristic algorithm explained in my ICC paper)
%
% -DEFINITIONS
%   1. 'bufferLinks' is same as "Broken Links" in the paper.
%
% -LOG
%   1. (2017/03/06): Changing 'requiredPowerScaling' from 1 to [1.01,1.1] improves the performance of 'sch_Heuristic1' from -2.5%  to +0.78% :-). However no performance difference for other schedulers.
%               requiredPowerScaling = 1.01 makes '#breakUponSatisfyingAllLinks' to break at 301th iteration, whereas requiredPowerScaling = 1.05 or 1.1 saves simulation time by breaking at 201th iteration itself
%
%
% -CONCLUSIONS
%   1. 2017/02/03: (20 VUE, single trial, natural-scheduling. Dummy-ACIR/LTE-ACIR): power-control improves #successfulLinks from 8 to 10. (ie. 25% improvement!!)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function powerPerTxPerPerTime = powerControl_kangWang(config,topology,txVUEPerRBPerTime,powerPerTxPerPerTime)


%% Parameter updation
config = patchConfig(config);
Ntx = config.N;
Nrx = config.N;
F = config.F;
T = config.T;
if(~any(txVUEPerRBPerTime))   
    powerPerTxPerPerTime = nan(Ntx,1,T);
    return;
end
if(gf.isempty('powerPerTxPerPerTime'))
    powerPerTxPerPerTime = 0.1*config.Pmax*matrixGeneralized(ones(Ntx,1,T), 0);
end
H = topology.H;
requiredPowerScaling = 1.05;

[A,b] = computeConstraintMatrixForPowerValues(topology,txVUEPerRBPerTime);
% A = gf.slack(A,2,1);
% b = b(:);

% -Definition:
%   link is basically for every timeslot. i.e. there are Ntx x Nrx x T links in total
% With 3D Matrix Assumption of  #RB  x  #RxVUEs  x  #Timeslots
% % % % lS.txVUEPerLink = gf.vec(repelem(txVUEPerRBPerTime,1,Nrx));
% % % % lS.rxVUEPerLink = gf.vec(repelem(1:Nrx,F,1,T));
% % % % lS.timeslotPerLink = gf.vec(repelem(shiftdim(1:T,-1),F,Nrx));
% -Trim-off non-scheduled links
lS = topology.getLinkStatistics(txVUEPerRBPerTime,A);
assert(all(lS.txVUEPerLink(lS.linksToConsider)>0) && all(lS.rxVUEPerLink(lS.linksToConsider)>0) ...
    && all(lS.timeslotPerLink(lS.linksToConsider)>0),'Error: lS.txVUEPerLink,lS.rxVUEPerLink,lS.timeslotPerLink must be nonzero integers!');
lC = lS.linksToConsider;
nLinks = numel(lS.txVUEPerLink);
deletionCountPerLink = zeros(nLinks,1);
requiredTxPowerPerLink = nan(nLinks,1);
for iter = 1:400
    deficientRxPowerPerLink = gf.vec(gf.mmat(A,powerPerTxPerPerTime) -b);
    requiredTxPowerPerLink(lC) = powerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),lS.txVUEPerLink(lC),1,lS.timeslotPerLink(lC))) + deficientRxPowerPerLink(lC) ./ H(sub2ind(size(H),lS.txVUEPerLink(lC),lS.rxVUEPerLink(lC)));
    requiredTxPowerPerLink = requiredPowerScaling*requiredTxPowerPerLink;       %% requiredPowerScaling 
    
    activeLinks = intersect(lC,find(deletionCountPerLink < config.kangWangConfig.Ct));
    % TERMINATION Condition
    successfulRealLinks = unique(lS.realLinkPerLink( activeLinks(deficientRxPowerPerLink(activeLinks)<=0) ));             %% Definition: 'realLink' corresponds to Ntx x Nrx links
    activeRealLinks = find(topology.findSuccessStatus(lS,activeLinks)==1);      
    if(isempty(setxor(successfulRealLinks,activeRealLinks)))                    %% #breakUponSatisfyingAllLinks
        break;
    end
    % -Debug
    % topology.findSuccessStatus(lS,activeLinks(deficientRxPowerPerLink(activeLinks)<=0))
    % Manage links
    bufferLinks = find(requiredTxPowerPerLink > config.Pmax);
    deletionCountPerLink(bufferLinks) = deletionCountPerLink(bufferLinks) +1;
    activeLinks = setdiff(activeLinks,bufferLinks);
    % Compute Tx-power
    % powerPerTxPerPerTime = computeTxPower_alg1(config,requiredTxPowerPerLink(activeLinks),lS.txVUEPerLink(activeLinks),lS.rxVUEPerLink(activeLinks),lS.timeslotPerLink(activeLinks));
    % powerPerTxPerPerTime = computeTxPower_alg2(config,requiredTxPowerPerLink(activeLinks),lS.txVUEPerLink(activeLinks),lS.rxVUEPerLink(activeLinks),lS.timeslotPerLink(activeLinks));
    powerPerTxPerPerTime = computeTxPower_alg3(config,requiredTxPowerPerLink(activeLinks),lS.txVUEPerLink(activeLinks),lS.rxVUEPerLink(activeLinks),lS.timeslotPerLink(activeLinks));
    assert(~gf.any(isnan(powerPerTxPerPerTime)),'Error: powerPerTxPerPerTime must not contain nan values');
end

end


%% Find optimal tx power
% -DEFINITIONS:
%       1. 'BLink'  ==  Bottlenecked Link
function powerPerTxPerPerTime = computeTxPower_alg1(config,requiredTxPowerPerLink,txVUEPerLink,rxVUEPerLink,timeslotPerLink)

    % -Parameter declarations
    Ntx = config.N;
    Nrx = config.N;
    F = config.F;
    T = config.T;
    powerPerTxPerPerTime = zeros(Ntx,1,T);
    requiredTxPowerPerTxPerRxPerTime = zeros(Ntx,Nrx,T);
    tempSize = size(requiredTxPowerPerTxPerRxPerTime);
    requiredTxPowerPerTxPerRxPerTime(sub2ind(tempSize,txVUEPerLink,rxVUEPerLink,timeslotPerLink))  =  requiredTxPowerPerLink;
    
    % -First allocate power for those tx-timeslot pairs having bottleneck links  #Algorithm
    [requiredTxPowerPerTxPerRxPerTime,requiredTxPowerPerTxPerPerTime] = updateRequiredPower(requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime);
    isBottleneckedPerTxPerRx = sum(requiredTxPowerPerTxPerRxPerTime>0,3)==1;      %% Bottleneck happens when a particular timesolot is mandatory for a Tx-Rx link  #Algorithm
    temp = nansum((requiredTxPowerPerTxPerRxPerTime>0) .* (bsxfun(@times,isBottleneckedPerTxPerRx,shiftdim(1:T,-1))),3);
    timeslotPerBLink = temp(temp>0);
    [txPerBLink,~] = find(isBottleneckedPerTxPerRx);
    if(~isempty(txPerBLink))
        powerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),txPerBLink,1,timeslotPerBLink)) = requiredTxPowerPerTxPerPerTime(gf.sub2ind(size(requiredTxPowerPerTxPerPerTime),txPerBLink,1,timeslotPerBLink));
    end

    % -Secondly find the tx-timeslot pairs which can serve maximum links, then allocate the power  #Algorithm
    for timeslotIter = 1:T
        [requiredTxPowerPerTxPerRxPerTime,requiredTxPowerPerTxPerPerTime] = updateRequiredPower(requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime);  % -Update requiredPower %% TODO: make this function
        nVUEsServedPerTxPerPerTime = sum(requiredTxPowerPerTxPerRxPerTime>0,2);         %% Number of (extra) VUEs served, if tx-VUE is transmitting on a timeslot
        [maxnVUEsServerPerTx,bestTimeslotPerTx] = max(nVUEsServedPerTxPerPerTime,[],3);
        bestTx = find(maxnVUEsServerPerTx);                                             %% Best transmitter to be allocated for current iteration
        if(isempty(bestTx)) break; end
        powerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),bestTx,1,bestTimeslotPerTx(bestTx))) = requiredTxPowerPerTxPerPerTime(gf.sub2ind(size(requiredTxPowerPerTxPerPerTime),bestTx,1,bestTimeslotPerTx(bestTx)));
    end
end




%% Find optimal tx power using the same algorithm presented in paper draft 2.2 of "Scheduling and Power control for V2V communication with ACI"
% -Definition:
%   "[requiredTxPowerPerTxPerRxPerTime]_{i,j,t} = nan"  =>  we don't need to consider the link (i,j,t)
function powerPerTxPerPerTime = computeTxPower_alg2(config,requiredTxPowerPerLink,txVUEPerLink,rxVUEPerLink,timeslotPerLink)

    % -Parameter declarations
    Ntx = config.N;
    Nrx = config.N;
    F = config.F;
    T = config.T;
    powerPerTxPerPerTime = zeros(Ntx,1,T);
    requiredTxPowerPerTxPerRxPerTime = nan(Ntx,Nrx,T);
    tempSize = size(requiredTxPowerPerTxPerRxPerTime);
    requiredTxPowerPerTxPerRxPerTime(sub2ind(tempSize,txVUEPerLink,rxVUEPerLink,timeslotPerLink))  =  requiredTxPowerPerLink;
    
    % -Algorithm-4 in paper draft 2.2
    while(~gf.all(isnan(requiredTxPowerPerTxPerRxPerTime)))
        [requiredTxPowerPerTxPerRx,bestTimeslotPerTxPerRx] = nanmin(requiredTxPowerPerTxPerRxPerTime,[],3);
        [requiredTxPowerPerTx,indPerTx] = nanmax(requiredTxPowerPerTxPerRx,[],2);
        tcapPerTx = bestTimeslotPerTxPerRx(gf.sub2ind(size(bestTimeslotPerTxPerRx),:,indPerTx));
        powerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),1:Ntx,1,tcapPerTx)) = requiredTxPowerPerTx;
        % Reset requiredTxPowerPerTxPerRxPerTime=nan for all those successful links
        [successfulTxs,successfulRxs,~] = gf.findIndices(bsxfun(@le,requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime));
        requiredTxPowerPerTxPerRxPerTime(gf.mysub2ind(size(requiredTxPowerPerTxPerRxPerTime),successfulTxs,successfulRxs,:)) = nan;
    end
    powerPerTxPerPerTime(isnan(powerPerTxPerPerTime)) = 0;
end


%% Find optimal tx power using the algorithm-4 presented in paper 'draft A4' of "Scheduling and Power control for V2V communication with ACI"
% -Definition:
%   "[requiredTxPowerPerTxPerRxPerTime]_{i,j,t} = nan"  =>  we don't need to consider the link (i,j,t)
function powerPerTxPerPerTime = computeTxPower_alg3(config,requiredTxPowerPerLink,txVUEPerLink,rxVUEPerLink,timeslotPerLink)

    % -Parameter declarations
    Ntx = config.N;
    Nrx = config.N;
    F = config.F;
    T = config.T;
    txVUEs = [1:Ntx]';
    powerPerTxPerPerTime = zeros(Ntx,1,T);
    requiredTxPowerPerTxPerRxPerTime = nan(Ntx,Nrx,T);
    tempSize = size(requiredTxPowerPerTxPerRxPerTime);
    requiredTxPowerPerTxPerRxPerTime(sub2ind(tempSize,txVUEPerLink,rxVUEPerLink,timeslotPerLink))  =  requiredTxPowerPerLink;

    for it = 1:T
        isSuccessPerTxPerRxPerTime = bsxfun(@le,requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime);
        requiredTxPowerPerTxPerRxPerTime(isSuccessPerTxPerRxPerTime) = nan;
        
        nLinksToServePerTxPerPerTime = sum(bsxfun(@ge,requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime),2);    %% Warning: Take care of NaN's while editing this line
        [nLinksToServeFortcapPerTx,tcapPerTx] = gf.maxRandIndex(nLinksToServePerTxPerPerTime,[],3);                 %% (2017 July 13): changed 'max' to 'gf.maxRandIndex', but hasn't verified the result change  #Log
        maxrequiredTxPowerPerTxPerPerTime = nanmax(requiredTxPowerPerTxPerRxPerTime,[],2);
        
        powerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),txVUEs(nLinksToServeFortcapPerTx>0),1,tcapPerTx(nLinksToServeFortcapPerTx>0))) ...
            = maxrequiredTxPowerPerTxPerPerTime(gf.sub2ind(size(powerPerTxPerPerTime),txVUEs(nLinksToServeFortcapPerTx>0),1,tcapPerTx(nLinksToServeFortcapPerTx>0)));
    end
end


%% If a link becomes successful with the existing 'powerPerTxPerPerTime', then flush (ie. zeroing) the corresponding value in 'requiredTxPowerPerTxPerRxPerTime'
function [requiredTxPowerPerTxPerRxPerTime,requiredTxPowerPerTxPerPerTime] = updateRequiredPower(requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime)
    isSuccessPerTxPerRxPerTime = bsxfun(@le,requiredTxPowerPerTxPerRxPerTime,powerPerTxPerPerTime);
    requiredTxPowerPerTxPerRxPerTime(isSuccessPerTxPerRxPerTime) = 0;
    requiredTxPowerPerTxPerPerTime = max(requiredTxPowerPerTxPerRxPerTime,[],2);
end



function config = patchConfig(config)
    %% -USER-DEFINED parameters:
    configDefault.probeLevel = 1;                   %% 1 => Capture Successful/Failed link indices, 2 => IoTPerHopPerLink Values
    configDefault.isDebug = false;
    configDefault.kangWangConfig.isDebug = false;
    %% -Create Default KangWang Config
    configDefault.kangWangConfig.removalPolicy = 'Distributed';            %% - {'Distributed','Centralised'}
    configDefault.kangWangConfig.termCondition = 'Feasiblity-Distributed'; %% - {'Convergence','Feasiblity','Feasiblity-Distributed'}
    configDefault.kangWangConfig.powerUpdationFn = @max;                   %% - {@max,@mean,@min}
    configDefault.kangWangConfig.linkRemovalScheme = 'kangWang';           %% - {'kangWang','bufferPercentage','bufferCount'}
    configDefault.kangWangConfig.Ct = 100;                                 %% - Buffer-Count threshold for 'kangWangModified' algorithm (Old variable 'bufferCountThreshold')
    %% -Merging Inputted config with Default inputs,
    config = gf.mergeTwoStructures(configDefault,config,'',true);
end



