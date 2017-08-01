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
%       1. Create a full graph with directiona? edges
%       2. Find the optimum vertex to remove, so that maximum outgoing edges are removed
%       3. Update the graph
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. 
%
%%%% -NOTES:
% 
%
%%%% -NOTES (Programming):
%
%%%% -TODO:
%       1.
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: varargin is placed so that this function won't throw error even while calling with more inputs 
function out = sch_edgeRemoval(config,~,~,~,blockInterleaverWidth,varargin)
    if(gf.isempty('blockInterleaverWidth'))
       blockInterleaverWidth = 3; 
    end
    N = config.N;
    F = config.F;
    T = config.T;
    edgePerTxVUEPerRxVUE = ones(N);
    edgePerTxVUEPerRxVUE(logical(eye(N))) = 0;
    %% First get the optimum list of VUEs to be scheduled in every time-slot 
    txVUEsPerTime = cell(1,T);
    for t = 1:T         %% For each Time-slot
        temp = edgePerTxVUEPerRxVUE;
        for f = 1:F     %% For each RB
            [~,scheduledVUE] = max(sum(temp,2));
            txVUEsPerTime{t}(end+1) = scheduledVUE;
            temp(scheduledVUE,:) = 0;
            temp(:,scheduledVUE) = 0;
        end
        edgePerTxVUEPerRxVUE(txVUEsPerTime{t},setdiff(1:N,txVUEsPerTime{t})) = 0;
    end
    %% Block interleave along freqency, in each time slot
    txVUEsPerTime = cellfun(@(x) sort(x),txVUEsPerTime,'UniformOutput',false);                                      %% Sort on every time-slot
    txVUEsPerTime = cellfun(@(x) gf.interleaveBlock_old(x,blockInterleaverWidth),txVUEsPerTime,'UniformOutput',false);  %% Do block interleaving
    out.txVUEPerRBPerTime = cell2mat(txVUEsPerTime);
end