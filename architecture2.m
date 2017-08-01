%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% -BRIEF DESCRIPTION:
%   
%
%%%% -DETAILED DESCRIPTION:
%       1. Output 'isSuccessPerTxPerRxPerAlgorithm' contains NaN, which indicates that either, no VUE transmitted on that frequency-time slot, 
%               or reception is blocked for transmitting VUEs when "config.isEnableFullDuplex = false;", or it's a self link (ie diagonal elements)
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. To compare nSuccessfulLinks across 'Algorithms'
%           tic, out1=architecture2('',15); toc
%           d = gf.congregate('out1(:).isSuccessPerTxPerRxPerAlgorithm');
%           gf.sum(d,[1,2,3])
%
%%%% -NOTES:
%       1. 'powerControl_kangWang' is same as 'Heuristic Power' in the paper
%       2. In both 'sch_naturalScheduling' & 'sch_blockInterleaver1',  a VUE is scheduled at max once.
%
%%%% -NOTES (Programming):
%       1. 'txVUE' = 0  =>  No VUE is transmitting
%       2. 'successStatus' = 1      =>  Success
%                          = 0      =>  Failure
%                          = NaN    =>  No attempted transmission or receptions
%       3. TERMINOLOGY UPDATION
%           a. 'sch_naturalScheduling'              <==>  BIS(w=1) in paper
%           b. 'sch_optimizedBlockInterleaver1'     <==>  BIS(optimized w) in paper
%           c. 'powerControl_kangWang'              <==>  Heuristic Power control in paper
%       4. Upon adding a new power-control algorithm, add exception to return "powerPerTxPerPerTime = nan(Ntx,1,T)" when no VUEs are scheduled.
%
%%%% -WARNINGS (Programming):
%       1. Matlab won't distinguish between a structure and, "an array of structers with numel '1' " 
%
%%%% -COMMANDS
%       1. tic, out1 = architecture2(); toc, squeeze(out1.txVUEPerRBPerTimePerAlgorithm), gf.nansum(out1.isSuccessPerTxPerRxPerAlgorithm,[1,2])
%       2. tic, outA = architecture2([],struct('nTrials',100)); toc, nanmean(gf.nansum(gf.congregate('outA(:).isSuccessPerTxPerRxPerAlgorithm')==1,[2,3]),1)
%
%
%%%% -OBSERVATIONS
%       1. (simulation of 2016 Dec 24): Order of scheduler performance:  'sch_naturalScheduling' < 'sch_blockInterleaver1' < 'sch_blockInterleaver3' < 'sch_blockInterleaver2' < 'sch_Heuristic1'  
%       2. 2016/02/03:  For single trial case, all scheduling algorithms work similar (ie. 8 & 10 successfulLinks without & with power control). The reason is that, the performance is limited by link-channel (not ACI).
%                           Corrected this by adding 'requiredPowerScaling' in 'powerControl_kangWang()'
%       3. 2016/02/09:  Doing scheduling & power-control again makes very little improvement in nSuccessfulLinks (6.3 to 6.4), but 14% reduction in total power!
%
%%%% -EXECUTION
%       1. 2016/02/03: All algorithms for single trial takes 1.5min
%       2. opt.schAlgorithmsToRun = {'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1','sch_Heuristic1('''','''','''','''',true)'};
%          'powerControl_equalPower' => , 'powerControl_gurobi' => 77 seconds
%
%
%%%% -TODO:
%       1.
%
%%%% -LOG
%       1. 2017/02/02: Changed Pmax from 30dBm to 24dBm.
%
%%%% -HashTags
%       1. '#Log' used to specify log
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function out = architecture2(config,opt)

% try 
% persistent architecture2_tryCatchDisplay; 
% if(isempty(architecture2_tryCatchDisplay)) architecture2_tryCatchDisplay = true;  disp('Warning: No try catch block in architecture2!!'); end

if(gf.isempty('config')) config = struct(); end
if(gf.isempty('opt.nTrials')) opt.nTrials = 1; end
if(gf.isempty('opt.schAlgorithmsToRun')) opt.schAlgorithmsToRun=''; end
if(gf.isempty('opt.powerControlAlgorithmsToRun')) opt.powerControlAlgorithmsToRun=''; end
if(gf.isempty('opt.isRecursiveCalling')) opt.isRecursiveCalling = false; end
if(~opt.isRecursiveCalling)
    % -Debug
%     rng('default'); disp('Warning: Resetting the seed');
    if(opt.nTrials==1) 
        config.configChannel.isShadow = false;
        config.isRandomDrop = false;
    else
        config.configChannel.isShadow = true;        
        config.isRandomDrop = true;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% - MANUAL INPUTS - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% config.T=1;
% config.N = 20;
% disp('T = 1  !!!');

if(gf.isempty('opt.schAlgorithmsToRun'))
    opt.schAlgorithmsToRun = {'sch_naturalScheduling'};%'sch_optimizedBlockInterleaver1','sch_Heuristic1'}; %%{'sch_naturalScheduling','sch_optimizedBlockInterleaver1','sch_Heuristic1','sch_gurobi'};   %%{'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1','sch_Heuristic1('''','''','''','''',true)'};   %%
%       opt.schAlgorithmsToRun = {'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1', ...  %% {'sch_naturalScheduling'};     %%  %% {'sch_gurobi'};       
%                         'sch_Heuristic1('''','''','''','''',2)','sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1','sch_gurobi'};
end

if(gf.isempty('opt.powerControlAlgorithmsToRun'))
    opt.powerControlAlgorithmsToRun =  {'powerControl_equalPower'}; %%{'powerControl_gurobi'};
%       opt.powerControlAlgorithmsToRun =  {'powerControl_equalPower','powerControl_equalPower','powerControl_equalPower','powerControl_equalPower','powerControl_equalPower', ...
%                        'powerControl_kangWang','powerControl_kangWang','powerControl_kangWang','powerControl_kangWang','powerControl_kangWang'};  %%{'powerControl_equalPower'};       %% 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


opt.powerControlAlgorithmsToRun = patchPowerControlAlgorithmsToRun(opt.powerControlAlgorithmsToRun,opt.schAlgorithmsToRun);
assert(numel(opt.powerControlAlgorithmsToRun) == numel(opt.schAlgorithmsToRun),'Error: #powerControlAlgorithmsToRun must be equals #schAlgorithmsToRun!!');

if(~opt.isRecursiveCalling)
    optCallee = opt;
    optCallee.isRecursiveCalling = true;
    optCallee.nTrials = 1;
    for i = 1:opt.nTrials
        out(i) = architecture2(config,optCallee);
        if(i==1)        %% Initializing out for speeding up simulation
            out(1:opt.nTrials,1) = out(1);
        end
    end
   return;
end


config = patchConfig(config);
Ntx = config.N;
Nrx = config.N;
F = config.F;
T = config.T;

if(F==0 || T==0 || Ntx==0 || Nrx==0)
    txVUEPerRBPerTimePerAlgorithm = [];
    powerPerTxPerTimePerAlgorithm = [];
    isSuccessPerTxPerRxPerAlgorithm = [];
    out = df.v2struct(txVUEPerRBPerTimePerAlgorithm,powerPerTxPerTimePerAlgorithm,isSuccessPerTxPerRxPerAlgorithm);
    return;
end

topology = Topology(config);
% H = topology.H;

txVUEPerRBPerTime_0 = [];
txVUEPerRBPerTimePerAlgorithm = nan(F,T,numel(opt.schAlgorithmsToRun));
powerPerTxPerTimePerAlgorithm = nan(Ntx,T,numel(opt.schAlgorithmsToRun));
isSuccessPerTxPerRxPerAlgorithm = zeros(Ntx,Nrx,numel(opt.schAlgorithmsToRun));
for iAlgorithm = 1:numel(opt.schAlgorithmsToRun)
    txVUEPerRBPerTime = [];
    powerPerTxPerPerTime = [];
    if(strcmp(opt.schAlgorithmsToRun{iAlgorithm},'sch_gurobi'))
        txVUEPerRBPerTime = txVUEPerRBPerTime_0;
    end
    schedulingCMD = gf.repairFunctionHandleString(opt.schAlgorithmsToRun{iAlgorithm},'config','topology','txVUEPerRBPerTime','powerPerTxPerPerTime');   %% Insert first 2 inputs to the schedulingCMD string
    for iter = 1:1
        % -Scheduling
        tempOut = eval(schedulingCMD);
        txVUEPerRBPerTime = tempOut.txVUEPerRBPerTime;
        % -Power Control
        powerPerTxPerPerTime = eval([opt.powerControlAlgorithmsToRun{iAlgorithm},'(config,topology,txVUEPerRBPerTime,powerPerTxPerPerTime)']);
%         %% Display Partial Results 
%         [~,isSuccessPerTxPerRxPerAlgorithm] = doTransmissionAndReception(topology,txVUEPerRBPerTime,powerPerTxPerPerTime);
%         disp(['iter = ',num2str(iter),'   Total number of successful links = ',num2str(gf.sum(isSuccessPerTxPerRxPerAlgorithm,[1,2])),'   mean power = ',num2str(gf.mean(powerPerTxPerPerTime,[1,2,3]))])
%         %% Nullify the scheduling incase zero power is allocated. (Useful only when iter>1)
%         [txVUEPerRBPerTime,powerPerTxPerPerTime] = refreshSchedulingAndPowerControlValues(txVUEPerRBPerTime,powerPerTxPerPerTime);        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -Result verification
    [~,isSuccessPerTxPerRxPerAlgorithm(:,:,iAlgorithm)] = topology.doTransmissionAndReception(txVUEPerRBPerTime,powerPerTxPerPerTime);
    % -Storing to output  (TODO_Future: Add scheduling and power control once more!)
    txVUEPerRBPerTimePerAlgorithm(:,:,iAlgorithm) = txVUEPerRBPerTime;
    powerPerTxPerTimePerAlgorithm(:,:,iAlgorithm) = permute(powerPerTxPerPerTime,[1,3,2]);
    if(strcmp(opt.schAlgorithmsToRun{iAlgorithm},'sch_Heuristic1'))        %% Saving the result for sch_gurobi
        txVUEPerRBPerTime_0 = txVUEPerRBPerTime;
    end
end

%% Make power values of non scheduled VUEs to NaN  (Not to 0, since nan is required in 'analyzeResult.majorAlgorithms')
powerPerTxPerTimePerAlgorithm = gf.preserveDataOnSpecifiedIndices(powerPerTxPerTimePerAlgorithm,txVUEPerRBPerTimePerAlgorithm,nan,1);
    
%% -Storing to output
out = df.v2struct(txVUEPerRBPerTimePerAlgorithm,powerPerTxPerTimePerAlgorithm,isSuccessPerTxPerRxPerAlgorithm);

%% -Debug
% display('#Successful Links: ');
% squeeze(gf.nansum(out.isSuccessPerTxPerRxPerAlgorithm,[1,2]))
% out.isSuccessPerTxPerRxPerAlgorithm(:,:,3)

% catch ME
%     ME = gf.saveWorkspaceUponException(ME);
%     throw(ME);
% end

end

%%
% -Definition:
%   VRB = Virtual frequency slots. 
function [out,txVUEPerVRBPerTime,freqSlotsCandidate] = sch_naturalScheduling(config,varargin)
    N = config.N;
    F = config.F;
    T = config.T;
    % -First Compute 'txVUEPerVRBPerTime' (ie. fill this F_ x T matrix)
%     if(~config.isEnableFullDuplex)      %% HALF-DUPLEX             (Comment(2017 June 30): The special case for fullDuplex is removed. this is for consistency while making report for the paper.)
%         if(T==1)
%             F_ = min(F,ceil(N/2));      % F_ is maximum number of frequency-slots, that can be scheduled in a timeslot
%         else
%             F_ = min(F,ceil(N/T));
%         end
%         candidateVUEs = gf.chooseMaximallySpreadElements([1:N]',F_*T);
%         txVUEPerVRBPerTime = gf.reshape(candidateVUEs,[T,F_],'-appendValue',0)';
%     else                                %% FULL DUPLEX
%        F_ = F;
%        candidateVUEs = gf.chooseMaximallySpreadElements([1:N]',F_*T);
%        txVUEPerVRBPerTime = reshape(gf.reshapeByCircularAffix(candidateVUEs,[F_*T,1]),[T,F_])';            
%     end
    
    if(~config.isEnableFullDuplex && T==1)
        F_ = min(F,ceil(N/2));      % F_ is maximum number of frequency-slots, that can be scheduled in a timeslot
    else
        F_ = min(F,ceil(N/T));
    end
    candidateVUEs = gf.chooseMaximallySpreadElements([1:N]',F_*T);
    txVUEPerVRBPerTime = gf.reshape(candidateVUEs,[T,F_],'-appendValue',0)';    
    % -Second step is to blow out matrix txVUEPerVRBPerTime (F_ x T)  to txVUEPerRBPerTime (F x T) matrix
%     w = (F-1)/(F_-1);
%     freqSlotsCandidate = round(1:w:F);
    freqSlotsCandidate = gf.chooseMaximallySpreadElements([1:F]',F_);
    out.txVUEPerRBPerTime = zeros(F,T);
    out.txVUEPerRBPerTime(freqSlotsCandidate,:) = txVUEPerVRBPerTime;
    out.txVUEPerRBPerTime = nullifyDuplicateElementsInAColumn(out.txVUEPerRBPerTime);   %% Added on 2016/12/12
    out.exitFlag = 1;
end

%% Same as the one used in architecture1
%   RBIndices [1:N] is filled in a block row by row, and read out column by column. 
%   'interleaverWidth': Adjacent VUEs would be using RBs 'interleaverWidth' apart  
function out = sch_blockInterleaver1(config,~,~,~,interleaverWidth)
    [out,txVUEPerVRBPerTime,freqSlotsCandidate] = sch_naturalScheduling(config);        % -Do natural Scheduling
    txVUEPerVRBPerTime = gf.interleaveBlock(txVUEPerVRBPerTime,interleaverWidth);       %% Interleave each columns
    out.txVUEPerRBPerTime(freqSlotsCandidate,:) = txVUEPerVRBPerTime;
    out.exitFlag = 1;
%     if(gf.isempty('interleaverWidth')) interleaverWidth = 3; end
%     if(~config.isEnableFullDuplex)
%         out = sch_blockInterleaver2(config,[],[],[],interleaverWidth,true);
%     else
%         N = config.N;
%         F = config.F;
%         T = config.T;
%         txVUEPerRBPerTime = rem(reshape(0:F*T-1,T,F)',N) +1;  
%         txVUEPerRBPerTime = nullifyDuplicateElementsInAColumn(txVUEPerRBPerTime);    
%         out.txVUEPerRBPerTime = gf.interleaveBlock(txVUEPerRBPerTime,interleaverWidth);  %% Interleave each columns
%         out.exitFlag = 1;
%     end
end

%%  Same as OBIS in the my paper "Scheduling and Power Control for Broadcast V2V Communications with Adjacent Channel Interference"
% Warning: Run 'run.findBestBlockInterleaverWidth1()' to create file 'Datafiles/sch_optimizedBlockInterleaver1.mat', before executing this function
function out = sch_optimizedBlockInterleaver1(config,~,~,~,interleaverWidth)
    persistent bestInterleaverWidthPerFsPerTsPerNs;
    if(isempty(bestInterleaverWidthPerFsPerTsPerNs))
        temp = load('Datafiles/sch_optimizedBlockInterleaver1.mat');
        bestInterleaverWidthPerFsPerTsPerNs = temp.bestInterleaverWidthPerFsPerTsPerNsPerDuplexTypePerACIRtype(:,:,:,config.isEnableFullDuplex+1,acir.acirIndexPerType(config.acirType));
    end

    interleaverWidth = bestInterleaverWidthPerFsPerTsPerNs(config.F,config.T,config.N);
    out = sch_blockInterleaver1(config,[],[],[],interleaverWidth);
end

%% - Commenting this function due to obsolescence, since this is now almost same as what Fredrik suggested.
% % %% Suggested by Fredrik on 2016/12/09
% % %   Fill 1:N VUEs row by row in 'txVUEPerRBPerTime', then interleave each column using standard interleaver of width 3. 
% % function out = sch_blockInterleaver2(config,~,~,~,interleaverWidth,isInterleaveAmongVUEs)        %%,interleaverWidth,varargin)
% %     if(gf.isempty('isInterleaveAmongVUEs')) isInterleaveAmongVUEs = true; end
% %     N = config.N;
% %     F = config.F;
% %     T = config.T;
% %     tmp = zeros(T,F);
% %     if(~config.isEnableFullDuplex && T==1)
% %         N_ = min(floor(N/2),F*T);       %% N_ is maximum number of VUEs, that can be scheduled in a timeslot
% %     else
% %         N_ = min(N,F*T);
% %     end
% %     schHopFactor = floor(N/N_);
% %     candidateVUEs = 1:schHopFactor:N;
% %     tmp(1:N_) = candidateVUEs(1:N_);
% %     txVUEPerRBPerTime = tmp';
% %     nTxVUEsInATimeSlot = sum(txVUEPerRBPerTime(:,1)>0);
% %     if(isInterleaveAmongVUEs)
% %         if(gf.isempty('interleaverWidth'))
% %             txVUEPerRBPerTime(1:nTxVUEsInATimeSlot,:) = gf.interleaveBlock(txVUEPerRBPerTime(1:nTxVUEsInATimeSlot,:),ceil(nTxVUEsInATimeSlot/2));   %% Interleave among VUEs
% %         else
% %             txVUEPerRBPerTime(1:nTxVUEsInATimeSlot,:) = gf.interleaveBlock(txVUEPerRBPerTime(1:nTxVUEsInATimeSlot,:),interleaverWidth);   %% Interleave among VUEs            
% %         end
% %     end
% %     out.txVUEPerRBPerTime = gf.interleaveBlock(txVUEPerRBPerTime,nTxVUEsInATimeSlot);                                                       %% Place 0s among Tx VUEs
% %     out.exitFlag = 1;
% % end


%% A VUE is trying to establish communication with 'interleaverWidth' on each side 
function out = sch_blockInterleaver3(config,~,~,~,interleaverWidth)        %%,interleaverWidth,varargin)
    error('TODO: Choose VUEs if N<FT ...');
    if(gf.isempty('interleaverWidth')) interleaverWidth = 3; end
    N = config.N;
    F = config.F;
    T = config.T;
    txVUEPerRBPerTime = rem(reshape(0:F*T-1,T,F)',N) +1;
    txVUEPerRBPerTime = nullifyDuplicateElementsInAColumn(txVUEPerRBPerTime);
    out.txVUEPerRBPerTime = gf.interleaveMaxSpread(txVUEPerRBPerTime);  %% Interleave each columns
    out.exitFlag = 1;
end


function powerPerTxPerPerTime = powerControl_equalPower(config,topology,txVUEPerRBPerTime,varargin)
    Ntx = config.N;
    T = config.T;
    if(~any(txVUEPerRBPerTime))   
        powerPerTxPerPerTime = nan(Ntx,1,T);
        return;
    end
    powerPerTxPerPerTime = config.Pmax*ones(Ntx,1,T);
end


%% Remove "zero powered VUEs" from 'txVUEPerRBPerTime'
%  Update "zero powered VUEs" power to reasonable values (This is for 2nd generation scheduling)
function [txVUEPerRBPerTime,powerPerTxPerPerTime] = refreshSchedulingAndPowerControlValues(txVUEPerRBPerTime,powerPerTxPerPerTime)
    [F,T] = size(txVUEPerRBPerTime);
    % -Refresh 'txVUEPerRBPerTime'
    powerPerRBPerTime = zeros(F,T);
    for t = 1:T
        [rbIndices,~,txVUEs] = find(txVUEPerRBPerTime(:,t));                %% Assumption 'txVUEPerRBPerTime' doesn't contain any NaN values
        powerPerRBPerTime(rbIndices,t) = powerPerTxPerPerTime(txVUEs,1,t);
    end
    txVUEPerRBPerTime = txVUEPerRBPerTime.*(powerPerRBPerTime>0);
    % -Refresh 'powerPerTxPerPerTime'
    powerPerTxPerPerTime_ = powerPerTxPerPerTime;
    powerPerTxPerPerTime_(powerPerTxPerPerTime_==0) = nan;
    avgPowerPerTime = shiftdim(nanmean(powerPerTxPerPerTime,1),2);
    for t = 1:T
        for iTx = find(isnan(powerPerTxPerPerTime_(:,1,t)))'
            temp = nanmean(powerPerTxPerPerTime_(iTx,1,:));     %% If iTx is transmitting in any timeslots, then take that power
            if(isnan(temp))                                     %% If iTx is not transmitting in no timeslots, then take average power across all transmitters in current timeslot
                temp = avgPowerPerTime(t);
            end
            powerPerTxPerPerTime_(iTx,1,t) = temp;
        end
    end
    powerPerTxPerPerTime = powerPerTxPerPerTime_;
end



% -Delete repeating elements in each column in 'txVUEPerRBPerTime' (The procedure is to take each row, and remove the repeating elements from below rows.) 
function txVUEPerRBPerTime = nullifyDuplicateElementsInAColumn(txVUEPerRBPerTime)
    F = size(txVUEPerRBPerTime,1);
    for f=1:F
        currentRow = txVUEPerRBPerTime(f,:);
        tmp = bsxfun(@minus,txVUEPerRBPerTime,currentRow);
        txVUEPerRBPerTime(tmp==0) = 0;
        txVUEPerRBPerTime(f,:) = currentRow;
    end
end

%%
%   If no 'powerControlAlgorithmsToRun' specified, then make it to 'powerControl_equalPower'
%   If 'powerControlAlgorithmsToRun' is singleton, then apply the same power control for all scheduling algorithms 
function powerControlAlgorithmsToRun = patchPowerControlAlgorithmsToRun(powerControlAlgorithmsToRun,schAlgorithmsToRun)
    if(gf.isempty('powerControlAlgorithmsToRun')) powerControlAlgorithmsToRun = cell(1,numel(schAlgorithmsToRun)); end
    if(numel(powerControlAlgorithmsToRun)==1) powerControlAlgorithmsToRun = repmat(powerControlAlgorithmsToRun,1,numel(schAlgorithmsToRun)); end
    for i=1:numel(powerControlAlgorithmsToRun)
        if(isempty(powerControlAlgorithmsToRun{i}))
            powerControlAlgorithmsToRun{i} = 'powerControl_equalPower';
        end
    end
end


%% -
function config = patchConfig(config)
    %% -Input Parameters
    configDefault.N = 20;    
    configDefault.F = 20;
    configDefault.T = 20;
    configDefault.isEnableFullDuplex = true;
    configDefault.configChannel.isNullifySelfChannel = true;
    configDefault.Pmax = db2pow(24-30);             %% Setting power to 24dbm
    configDefault.gammaT = db2pow(5);
    configDefault.noiseVariance = 12*db2pow(-136);  %% Noise Power for single RB
    configDefault.configChannel.isShadow = false;   %% TOOD
    configDefault.acirType = 'LTEmask';     %% 'LTEmask','scFDMA'
    configDefault.isRandomDrop = false;
    configDefault.adjacentVehicularDistance = 48.6;     %% 48.6m corresponds to 2.5sec for 70Km/hr as specified in 3GPP 36.885 section A.1.2
    configDefault.adjacentVehicularDistanceMin = 10;     %% (Obsolete note:We set min distance to 0, and lowerbound all distance to 10m while computing pathloss using Karedal model)
    %% -Merging Inputted config with Default inputs,
    config = gf.mergeTwoStructures(configDefault,config,'',true);
end
