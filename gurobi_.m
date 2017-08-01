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
%       Same as gurobi function, by providing some default inputs
%
%%%% -DETAILED DESCRIPTION:
%       1.
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



function [x,exitFlag,result] = gurobi_(model,params)

    [model,params] = patchConfig('model','params');

    result = gurobi(model,params);
    x = result.x;
    exitFlag = mapGurobiStatusToexitFlag(result.status);

end

%% Convert Gurobi status string to  Matlab linprog exitFlag 
%   Gurobi Status messages: http://www.gurobi.com/documentation/6.0/refman/optimization_status_codes.html#sec:StatusCodes  
%   Matlab exitFlags      : https://se.mathworks.com/help/optim/ug/linprog.html#outputarg_exitflag  
function exitFlag = mapGurobiStatusToexitFlag(gurobiStatus)
    switch(gurobiStatus)
        case 'OPTIMAL'
            exitFlag = 1;
        case 'INFEASIBLE'
            exitFlag = -2;
        case {'ITERATION_LIMIT','TIME_LIMIT'}
            exitFlag = 0;
        otherwise
            exitFlag = 100;
    end
end
    
    
%% -
function [model,params] = patchConfig(modelString,paramsString)
    %% Bring inputs from caller's scope
    if(evalin('caller',['gf.isempty(''',modelString,''')'])) model = struct();
    else model = evalin('caller',modelString); end;
    if(evalin('caller',['gf.isempty(''',paramsString,''')'])) params = struct();
    else params = evalin('caller',paramsString); end;
    %% -Default Parameters
    modelDefault = struct();  paramsDefault = struct();
    modelDefault.sense = '<';
    modelDefault.modelsense = 'min';
% % %     modelDefault.start = x0;
% % %     paramsDefault.TimeLimit = 300;
%     paramsDefault.method = 2;              %% Barrier method REF: http://www.gurobi.com/documentation/6.0/refman/method.html#parameter:Method
%     paramsDefault.outputflag = 0;
%     paramsDefault.resultfile = 'gurobi.log';
%     paramsDefault.Cuts = 3;
%     paramsDefault.MIPFocus = 1;            %% Prioritise getting better feasible solution, rather than reducing optimality-gap
%     paramsDefault.ImproveStartTime = 30;   %% Focus on getting better feasible solution after this time
%     paramsDefault.Heuristics = 1;
%     paramsDefault.MIPGap = 0.10;           %% 10% relative gap
    %% -Merging Inputted config with Default inputs,
    model = gf.mergeTwoStructures(modelDefault,model,'',true);
    params = gf.mergeTwoStructures(paramsDefault,params,'',true);
end    
    
    
    
    
    
    
