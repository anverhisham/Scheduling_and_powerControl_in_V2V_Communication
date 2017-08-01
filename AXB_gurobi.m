%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%  D = | kron(B',A),                I(N^2,N^2)  |
%      | col& row sum inequality,   0(2N,N^2)   | 
%      | total sum inequality,      0(1,N^2)    |  
%
%  TOOD:
%   Getting following error message:
%       Error using gurobi_mex
%       Error: could not create environment. Please make sure that GUROBI has a license and been validated.
%
%  Obervation:
%   1. Only 3% chance for getting 6th bit error (real outcome is better than the predicted by Gurobi)
%       Remaining errors are having significantly less chance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function out = AXB_gurobi(A,B,gamma)

    out.exitStatus = 0;
    t1 = cputime;
    F = size(A,2);      %% Number of RBs
    N = size(B,1);
    
    if(N>15)
        out.X = eye(N);
        out.simulationTime = nan;
        out.exitStatus = 2^10;
        return;
    end
    
    
    %% -Convert AXB problem into min c'x s.t. Dx <= b problem, without boolean constraints
    [c,D,b] = AXBproblem_to_Dxb_translate(A,B,gamma);
    nVars = numel(c);
    
    %% -Create a feasible starting point xAndy0  (Conclusion: Following xAndy0 increases simulation time by ~7 times!)
    X0 = eye(N);    
    x0 = X0(:);
    y0 = gf.vec(A*X0*B < gamma);
    xAndy0 = [x0;y0];
% %     %% -Standard Gurobi setups (Borrowed from test_gurobi_mex_MIP.m)
% %     opts.Start = xAndy0;
% %     vtypes = repmat('B',1,nVars);
% %     lb = zeros(nVars,1); ub = ones(nVars,1);
% %     contypes = '<';
% %     objtype = 1; % 1 - minimize, -1 - maximize
% %     clear opts;
% %     opts.IterationLimit = 1e7;
% %     opts.FeasibilityTol = 1e-9;
% %     opts.IntFeasTol = 1e-9;
% %     opts.OptimalityTol = 1e-9;
% %     opts.Method = 1;    % 0 - primal, 1 - dual
% %     opts.Presolve = -1; % '-1' - auto, 0 - no, 1 - conserv, 2 - aggressive
% %     opts.Display = 0; opts.DisplayInterval=0; opts.OutputFlag=0; opts.Display=0;
% %     opts.LogFile = 'AXB_gurobi.log';
% %     opts.WriteToFile = 'AXB_gurobi.mps';
% %     %  -Running Optimization tool
% %     [xAndy,val,exitflag,output] = gurobi_mex(c,objtype,sparse(D),b,contypes,[],[],vtypes,opts);    
    
    %% -Another way of setting parameters (Ref: http://www.gurobi.com/documentation/6.0/quickstart_mac/matlab_example.html)
    clear model;
    model.A = sparse(D);
    model.obj = c;
    model.rhs = b;
    model.sense = '<';
    model.vtype = 'B';
    model.modelsense = 'min';
    model.start = xAndy0;
    clear params;
    params.TimeLimit = 300;
    params.method = 2;              %% Barrier method REF: http://www.gurobi.com/documentation/6.0/refman/method.html#parameter:Method
    params.outputflag = 0;
    params.resultfile = 'AXB_gurobi.log';
    params.Cuts = 3;
    params.MIPFocus = 1;            %% Prioritise getting better feasible solution, rather than reducing optimality-gap
    params.ImproveStartTime = 30;   %% Focus on getting better feasible solution after this time
    params.Heuristics = 1;
    params.MIPGap = 0.10;           %% 10% relative gap
%     params.TimeLimit = 600;   %% -Specify max time-limit

    result = gurobi(model, params);
    xAndy = result.x;
    
    %% -Now Verify the result
    if(strcmpi(result.status,'OPTIMAL'))
        x = xAndy(1:N^2); y = xAndy(N^2+1:end);
        C1 = D(1:N^2,1:N^2);    %% -Same as C1 in function 'AXBproblem_to_Dxb_translate'
        if(sum(C1*x<=b(1:N^2))>sum(y==0))           %% 6th bit is set if the real outcome is better than the predicted by Gurobi
            out.exitStatus = 2^5;
        elseif(sum(C1*x<=b(1:N^2))<sum(y==0))       %% 7th bit is set if the real outcome is worse than the predicted by Gurobi
            out.exitStatus = 2^6;
        end
    else
       x = eye(N);
       if(strcmpi(result.status,'ITERATION_LIMIT') || strcmpi(result.status,'NODE_LIMIT') || strcmpi(result.status,'TIME_LIMIT'))
           out.exitStatus = 2^10;
       else
           out.exitStatus = 2^14;
           out.solverSay = result;
       end
    end

    %% -
    out.X = reshape(x,N,N);
    out.simulationTime = cputime - t1;
    
end
