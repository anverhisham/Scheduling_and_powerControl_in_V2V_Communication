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
%   N -> #VUEs
%   F -> #RBs
%   T -> #Time-slots
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. 
%
%%%% -NOTES:
%       1. Full Duplex reception is disabled/enabled using 'isSuccessPerRBPerRxVUE(:,txVUEPerRB) = 0;'
%
%%%% -NOTES (Programming):
%
%%%% -TODO:
%       1. Make listening-only VUEs to include in 'isSuccessPerRBPerTimePerRx' 
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Topology

properties
    config;
    locationPerNode;
    H;          %% NxN Channel matrix
    p;          %% Nx1 vector
    Lambda;     %% FxF ACI matrix
    I;          %% Generalised NxN identity matrix 
end

methods

    %% Inputs
    %   N - number of VUEs
    function this = Topology(config)
        this.config = config;
        this.locationPerNode = dropUEs(this,config.N,config);
        this.H = computeALLtoALLChannelGain(this.locationPerNode,config.configChannel);
        this.Lambda = acir.getACIRfromInterfererToDesiredLink([1:config.F],config.acirType);
        this.I = matrixGeneralized(eye(config.N),0);
        this.p = config.Pmax*ones(config.N,1);      %% TODO_Future: make separate function for setting power values.
    end

    %% Notes:
    %       1. while using input 'txVUEPerRBPerTime':
    %               Use '0' for non-scheduled fequency-time slot in input 'txVUEPerRBPerTime' 
    %       2. txVUEPerRBPerTime = All zero matrix   =>  isSuccessPerRBPerTimePerRx = all NaN matrix
    function [isSuccessPerRBPerTimePerRx,isSuccessPerTxPerRx] = doTransmissionAndReception(this,txVUEPerRBPerTime,powerPerTxPerPerTime)
        Ntx = this.config.N;
        Nrx = this.config.N;
        F = size(txVUEPerRBPerTime,1);              %% F is the number of RBs
        T = size(txVUEPerRBPerTime,2);        
        if(gf.isempty('powerPerTxPerPerTime')) powerPerTxPerPerTime = repmat(this.p,1,1,T); end        
        [A,b] = computeConstraintMatrixForPowerValues(this,txVUEPerRBPerTime);
        isSuccessPerRBandRxVUEPerPerTime = double(gf.mmat(A,powerPerTxPerPerTime) <= b);
        % -Rearrange to make to final output
        temp = reshape(isSuccessPerRBandRxVUEPerPerTime,F,Nrx,T);
        isSuccessPerRBPerTimePerRx = permute(temp,[1,3,2]);
        % -If no VUE is scheduled on a time-frequency slot, then make success status nan,
        isSuccessPerRBPerTimePerRx(repmat(~(txVUEPerRBPerTime>0),[1,1,Nrx])) = nan;     %% '~' is used to include 'nan' also
        if(~this.config.isEnableFullDuplex)                                             %% Ignore jth VUE reception link (by making it NaN), if jth VUE is already scheduled in any of the RBs.
            txVUEVectorPerRBPerTime = repmat(permute(txVUEPerRBPerTime,[3,2,1]),[F,1]);
            isSuccessPerRBPerTimePerRx = gf.clearDataOnSpecifiedIndices(isSuccessPerRBPerTimePerRx,txVUEVectorPerRBPerTime,nan,3);
        end
        % - Remove success status of self-channel
        [rbIndices,timeslots] = find(txVUEPerRBPerTime);
        txVUEs = txVUEPerRBPerTime(txVUEPerRBPerTime>0);
        isSuccessPerRBPerTimePerRx(sub2ind(size(isSuccessPerRBPerTimePerRx),rbIndices,timeslots,txVUEs)) = nan;
        if(nargout>=2)
           isSuccessPerTxPerRx = sf.getTxRxSuccessStatusFromReceptionStatus(txVUEPerRBPerTime,isSuccessPerRBPerTimePerRx); 
        end
    end
    
    %%
    %   Definition:
    %       A = (Nrx*F x Ntx x nTimeslots) matrix (with first F rows coressponds to 1st VUE reception. For example, (F+3)rd row corresponds to 2nd VUE reception on 3rd RB)
    %   Note:
    %       1. 3rd dimension of A & b is for timeslot
    %   Convert    X*P*H - gammaT*Lambda*X*P*H >= gammaT*noiseVariance     to     Ap  <=  b
    %   Look at fair-notebook (2017/2/2) for more understanding
    function [A,b] = computeConstraintMatrixForPowerValues(this,txVUEPerRBPerTime)
        Lambda = this.Lambda;
        Ntx = this.config.N;
        Nrx = this.config.N;
        F = size(txVUEPerRBPerTime,1);
        T = size(txVUEPerRBPerTime,2);
        gammaT = this.config.gammaT;
        noiseVariance = this.config.noiseVariance;
        I = eye(F);
        I_Ntx = eye(Ntx);
        H = this.H;
        Hexpanded = repelem(H.',F,1); 
        A = nan(Nrx*F,Ntx,T);  b = nan(Nrx*F,1,T);
        for t = 1:T                                 %% For each time-slot
            txVUEPerRB = txVUEPerRBPerTime(:,t);
            % -Create F x Ntx Scheduling matrix 'X'
            % X = this.I(txVUEPerRB,:);                    %% X_{f,i} indicate if fth RB is scheduled to ith VUE  %% (2017/02/23): Remove this for speed
            [RBs,~,txVUEs] = find(txVUEPerRB);
            X = zeros(F,Ntx);
            X(RBs,:) = I_Ntx(txVUEs,:);
            A(:,:,t) = -repmat((I-gammaT*Lambda)*X,Nrx,1).*Hexpanded;   %% lth row of 'A' is the SINR constraint for the VUE receiver "floor(l,Nrx)" on RB "l%F"
        end
        b(:) = -gammaT*noiseVariance;
        % -Remove self links, by making the corresponding contraint impossible
        rxVUEPerLinkPerPerTime = repelem([1:Nrx]',F,1,T);         %% There are F x Nrx links in a timeslot
        txVUEPerLinkPerPerTime = repmat(permute(txVUEPerRBPerTime,[1,3,2]),Nrx,1);
        for t = 1:T
            A(rxVUEPerLinkPerPerTime(:,:,t)==txVUEPerLinkPerPerTime(:,:,t),:,t) = Inf;
        end
    end
    
    
    %% - Input 'N' is same as number of UEs ...
    function locationPerUE = dropUEs(this,N,config)
        if(~config.isRandomDrop)
            locationPerUE = 0:config.adjacentVehicularDistance:(N-1)*config.adjacentVehicularDistance;
        else
            locationPerUE = cumsum([0;exprnd(config.adjacentVehicularDistance-config.adjacentVehicularDistanceMin,N-1,1)+config.adjacentVehicularDistanceMin]);
        end
    end

    
    function linkStatistics = getLinkStatistics(this,txVUEPerRBPerTime,A)
        config = this.config;
        Ntx = config.N;
        Nrx = config.N;
        F = config.F;
        T = config.T;
        % -Get all parameters related to currentLinks
        [rbPerLink,rxVUEPerLink,timeslotPerLink] = gf.findIndices(ones(F,Nrx,T));
        txVUEPerLink = gf.vec(repmat(permute(txVUEPerRBPerTime,[1,3,2]),1,Nrx));
        realLinkPerLink = (rxVUEPerLink-1)*Ntx + txVUEPerLink;      %% Definition: 'realLink' corresponds to Ntx x Nrx links
        linksToConsider = find(gf.vec(any(A<0,2)));                 % -Trim-off non-scheduled currentLinks
        % -Ignore reception links for transmitters, if full-duplex
        if(config.isEnableFullDuplex==false)
            for t=1:T
                links = find(timeslotPerLink==t);
                txVUEs = unique(txVUEPerLink(links));
                linksViolatingHalfDuplex = links(ismember(rxVUEPerLink(links),txVUEs));
                linksToConsider = setdiff(linksToConsider,linksViolatingHalfDuplex);
            end
        end
        linkStatistics = df.v2struct(rbPerLink,rxVUEPerLink,timeslotPerLink,txVUEPerLink,realLinkPerLink,linksToConsider);
    end
    
    
    function isSuccessPerTxPerRx = findSuccessStatus(this,linkStatistics,isSuccessPerLinkOrSuccessfulLinks)
        config = this.config;
        Ntx = config.N;
        Nrx = config.N;        
        F = config.F;
        T = config.T;
        if(islogical(isSuccessPerLinkOrSuccessfulLinks))
            isSuccessPerLink = isSuccessPerLinkOrSuccessfulLinks;
        else
            isSuccessPerLink = false(F*Nrx*T,1);
            isSuccessPerLink(isSuccessPerLinkOrSuccessfulLinks) = 1;
        end
        isSuccessPerTxPerRx = nan(Ntx,Nrx);
        isSuccessPerLink = logical(isSuccessPerLink);
        isSuccessPerTxPerRx(linkStatistics.realLinkPerLink(isSuccessPerLink)) = 1; 
    end
    
    function links = convertRealLinksToLinks(this,linkStatistics,realLinks)
        links = find(ismember(linkStatistics.realLinkPerLink,realLinks));
    end
    
    
        
end
end