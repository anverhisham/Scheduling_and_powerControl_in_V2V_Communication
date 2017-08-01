%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Co
    
properties (Constant)
    %% -Global Rules:
    indexPerSchAlgorithm = containers.Map({'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1','sch_Heuristic1([],[],[],[],true)'...
            'sch_gurobi','sch_gurobi('' '','' '',txVUEPerRBPerTime_0)','sch_gurobi(config,topology,txVUEPerRBPerTime_0)','sch_optimizedBlockInterleaver1'},{1,2,3,4,5,6,7,7,7,8});  %% NOTE: 'sch_optimizedBlockInterleaver1' is same as OBIS in the paper
    indexPerPowerControlAlgorithm = containers.Map({'powerControl_equalPower','powerControl_kangWang','powerControl_gurobi'},{1,2,3});      %% In case, changing this order, look for #AffectedByPCIndex
     
end

end
