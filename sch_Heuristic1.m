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
%       1. For each time-frequency slot
%           1.1 Compute nSuccessfulLinks upon scheduling each VUE   
%           1.2 Choose the best VUE, and schedule it.
%       2. 'iteration'  =>  Number of sweeping done. One can control the path of sweep using 'schedulingOrderPerGrid'
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. ans1=[probe.sch_Heuristic1.nSuccess]; ans2=reshape(ans1,2,[]); figure; plot(ans2(2,:)-ans2(1,:)); mean(ans2(2,:)-ans2(1,:))
%
%%%% -NOTES:
% 
%
%%%% -NOTES (Programming):
%       1. VUE=0  =>  No VUE is transmitting!
%
%%%% -STUDY ITEMS:
%       1. Find out how many iteration required for convergence.
%       
%
%%%% -OBSERVATION (Without Shadowing):
%       1. rbSchedulingSequence=gf.interleaveMaxSpread([1:F]')  reduces performance. Hence use natural order (i.e. [1,2,3,...])
%       2. Should scheduling done column by column (ie, first along RBs) or Row by Row ?
%               Column by Column scheduling seems to provide better performance
%       3. Do iter have any effect?
%               No effect, even when we change the order of scheduling during 2nd iteration.
%       4. (2016/12/11): 2nd iteration Row sweeping's result is stored today. row-sweeping iteration improves the result (#successfulLinks) by ~0.2.
%           (Improvement is highest when F=20 & T=2)
%       5. (2016/12/11): 2nd iteration Column sweeping's result is stored today. col-sweeping result is same as above row-sweeping
%       6. Let's try 3rd iterative sweeping
%
%%%% -TODO:
%       1. Create a Universal probe class, to probe why no improvement over 2nd iteration.
%       2. For tie-breaking, use surplus SINR or something...
%       3. Let Sweep start from middle RB, instead of end RB
%       4. Schedule 2RBs at a time, instead of 1RB
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = sch_Heuristic1(config,topology,txVUEPerRBPerTime,powerPerTxPerPerTime,algorithmIndex)

    if(gf.isempty('algorithmIndex')) algorithmIndex = 1; end
    persistent Fold Told sweepOrder_1 sweepOrder_2;

    global probe;
    Ntx = config.N;
    Nrx = config.N;
    F = config.F;
    T = config.T;
    tol = 1e-13;        %% NOTE: Since quantization error depends upon the order of addition matters, we need this tolerance while comparing numbers
     
    isSuccessPerRBPerTimePerRx = nan(F,T,Nrx);
    if(gf.isempty('powerPerTxPerPerTime')) powerPerTxPerPerTime = config.Pmax*ones(Ntx,1,config.T); end
    if(gf.isempty('txVUEPerRBPerTime')) txVUEPerRBPerTime = zeros(F,T); end
    nSuccess = 0;
    
    vueTuples = num2cell([0:Ntx]');
    if(algorithmIndex==1)
        %sweepOrder_1 = reshape([1:F*T],F,T);  
        %sweepOrder_2 = reshape([1:F*T],T,F)'; 
        if(Told==T & Fold==F)
            % Do nothing
        else
            fh = @(x) -sf.findACIInAnRB(x,config.acirType);
            sweepOrder_1 = sf.findSchedulingRBSequence(F,fh);                       %%  ** Here we are finding the scheduling order f_arrow ** %%
            sweepOrder_1 = bsxfun(@plus,sweepOrder_1(:,1),[0:T-1]/(2*F*T));
            sweepOrder_2 = bsxfun(@plus,sweepOrder_1(:,1),[0:T-1]*F);
            Fold = F;   Told = T;
        end
    elseif(algorithmIndex==2)
        vueTuples = num2cell(nchoosek([0:Ntx],2),2);
        vueTuples = [vueTuples; cellfun(@fliplr,vueTuples,'UniformOutput',false)];
        vueTuples = [{[0,0]};vueTuples];
        sweepOrder_1 = [9,9,7,7,5,5,3,3,1,1,2,2,4,4,6,6,8,8,10,10]'; sweepOrder_1 = bsxfun(@plus,sweepOrder_1,[1:F:F*T]);
        sweepOrder_2 = bsxfun(@plus,T*sweepOrder_1(:,1),[1:T]);
    elseif(algorithmIndex==3)       %% Instead of starting from end RB, start from middle RB
        sweepOrder_1 = [19,17,15,13,11,9,7,5,3,1,2,4,6,8,10,12,14,16,18,20]';
        sweepOrder_1 = bsxfun(@plus,sweepOrder_1,[0:T-1]*F);
        sweepOrder_2 = bsxfun(@plus,sweepOrder_1(:,1),[0:T-1]/(2*F*T));
    end
    assert(all([F,T]<=size(sweepOrder_1)),'Error: Sweep order not defined for high values of F & T !!');
    schedulingOrderPerGrid{1} = sweepOrder_1(ceil(end/2) - ceil(F/2)+1:ceil(end/2) + floor(F/2), :);
    schedulingOrderPerGrid{2} = sweepOrder_2(ceil(end/2) - ceil(F/2)+1:ceil(end/2) + floor(F/2), :);
    schedulingOrderPerGrid{3} = schedulingOrderPerGrid{1}; 
    
    linkWeightMatrix = getLinkWeightMatrix(ones(Ntx,Nrx));
    for iter = 1:1      %% Observation (2017/03/13): 3 iterations (compared to 1) improves the result by 0.18% only, but increases runtime from 1.3s to 3.9s for [F,T,N] = [20,5,20] for 100trials
        for ind = sort(unique(schedulingOrderPerGrid{iter}(schedulingOrderPerGrid{iter}>0)))'
            gridIndices = find(schedulingOrderPerGrid{iter}==ind);
            [fs,ts] = find(schedulingOrderPerGrid{iter}==ind);
            isSuccessPerRBPerTimePerRx_ = isSuccessPerRBPerTimePerRx;
            txVUEPerRBPerTime_ = txVUEPerRBPerTime;
            for n=1:numel(vueTuples)
%                 if(n~=0 && any(txVUEPerRBPerTime(:,ts)==n)) %% Avoid repeated scheduling of a VUE in a timeslot
%                    continue;
%                 end
                vueTuples_ = gf.reshape(vueTuples{n},size(txVUEPerRBPerTime(gridIndices)));   %% if(numel(gridIndices)<numel(vueTuples{n})), then trimm off 'vueTuples{n}'
                if(gf.all(txVUEPerRBPerTime(gridIndices)==vueTuples_))          %% Avoid scheduling the already scheduled VUEs in these RB-timeslots
                   continue;
                end
                txVUEPerRBPerTime_(gridIndices) = vueTuples_;
                isSuccessPerRBPerTimePerRx_(:,ts,:) = topology.doTransmissionAndReception(txVUEPerRBPerTime_(:,ts),powerPerTxPerPerTime(:,:,ts));
                isSuccessPerTxPerRx = sf.getTxRxSuccessStatusFromReceptionStatus(txVUEPerRBPerTime_,isSuccessPerRBPerTimePerRx_);
                %nSuccess_ = gf.sum(isSuccessPerTxPerRx==1,[1,2]);
                nSuccess_ = nansum(isSuccessPerTxPerRx(:).*linkWeightMatrix(:));        %% 'linkWeightMatrix' is brought in, instead of simply counting the number of successful links, to break the tie.
                if(nSuccess_ > nSuccess +tol)     %% Observation: (2017/03/12: Adding randomness here worsen the performance for single or 100 trials) 
                    nSuccess = nSuccess_;
                    isSuccessPerRBPerTimePerRx = isSuccessPerRBPerTimePerRx_;
                    txVUEPerRBPerTime = txVUEPerRBPerTime_;
                end
                isSuccessPerTxPerRx_bk = isSuccessPerTxPerRx;
            end
        end
% %         %%%% PROBING  (NOTE(2017/02/23): Commented it, since gf.isempty is taking 275s when algorithmIndex=true !!)
% %         if(~gf.isempty('probe.isDebug') && probe.isDebug)
% %             if(gf.isempty('probe.sch_Heuristic1'))
% %                 probe.sch_Heuristic1(1) = struct('Ntx',Ntx,'F',F,'T',T,'nSuccess',nSuccess);
% %             else
% %                 probe.sch_Heuristic1(end+1) = struct('Ntx',Ntx,'F',F,'T',T,'nSuccess',nSuccess);
% %             end
% %         end
    end
    out.txVUEPerRBPerTime = txVUEPerRBPerTime;
    out.exitFlag = 1;
end

%% Make NxN weight matrix, preferring neighbouring links over far away links
function linkWeightMatrix = getLinkWeightMatrix(currentWeightMatrix)
    [Ntx,Nrx] = size(currentWeightMatrix);
    scalingFactor = 1/(numel(currentWeightMatrix));
    linkWeightMatrix = currentWeightMatrix + scalingFactor*toeplitz(2.^-[1:Ntx],2.^-[1:Nrx]);
end








