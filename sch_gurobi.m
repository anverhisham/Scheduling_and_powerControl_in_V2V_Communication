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
%       See the fair-notebook for detailed description
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       
%
%%%% -NOTES:
%       1. X = FxNtx matrix with "X(f,i) = 1", if ith UE is scheduled on fth RB
%          AXB = FxNrx matrix with "AXB(f,j) >= gamma", if packet transmitted on fth RB is received successfully on jth receiver
%       2. It looks like halfDuplex not possible, by allowing self channel 'H_ii' a high values, since
%           a). There is high dynamic range b/w A & b. Each row of A would contain very high values (both self channel & eta) ~10^6, while b~=1.
%           b). To avoid the above problem, one approach is to reduce the self channel, and minimize eta by using "A = upperboundHighestValuesInA(A,b)".
%                But this has got following problem; reducing H_ii would nullify "aci*H_ii", hence allowing full duplex when tx & rx RBs are far apart.
%          To overcome this problem, I have transformed problem from "AXB + eta*Y >= gamma" to "AXB + 2*eta*Y - eta*ones(F,F)*X >= gamma"
%                   
%
%%%% -NOTES (Programming):
%       
%
%%%% -STUDY ITEMS:
%       1. TOOD_Urgent: [F,N,T]=deal(4,3,1); X = sym('X',[F,N,T]); H = sym('H',[N,N,T]); A_=-10.^-[0,1,2,3; 1,0,1,2; 2,1,0,1; 3,2,1,0]; A_(logical(eye(size(A_)))) = 1; symOut=A_*X*H;  why symOut(1,2) looks strange??
%       
%
%%%% -OBSERVATION (Without Shadowing):
%       1. Simulation is NOT POSSIBLE for T=3, due to high simulation time!! When T=3, there 28000 symbolic variables, ~20000 inequalities. There are 2 bottlenecks
%           a. Concatenating large symbolic matrices takes too long time. (However this can be avoided by computing and concatenating coefficients instead)
%           b. Finding coefficients take enormous time (~2min per inequality)! Inside 'gf.findCoefficientsPerVariable', the bottleneck is finding the location of observed coefficients in the list of all 28000 coefficients using ismember()
%         Possible Way Forwards;
%           i. Try to redesign the system with less number of symbolic variables, and constraints
%           ii. Optimize the system only for the case T=1, by avoiding variables constraining across timeslots!
%       2. (2017/02/14): [N,F,T] = [20,20,1],  =>  #successfulLinks for ['alg_Natural','sch_Heuristic1'] = [74,100]
%           params.Method = Default
%               Taking long running time; Optimality gap is (5.1% without starting point, 5.08% with starting point) even after 2hours for 
%               after that it's taking 10min for every 0.01% improvement! :-(.  Update after 4.52%, it's taking hours for each 0.01% increment!!
%               After 30seconds simulation [#successfulLinks,UpperBound] = [105,138], after 10min simulation [105,124.3], and after 15 hours it's [106,119.25]  =>  No improvement with hours of simulation!
%           params.Method = 1 ('Dual-simplex')
%               After 30seconds simulation [#successfulLinks,UpperBound] = [105,138], after 5min [105,124.5]
%           params.Method = 2 ('Barrier')
%               After 30seconds simulation [#successfulLinks,UpperBound] = [106,137.9], after 5min [106,123.9]
%           params.Method = 3 ('Concurrent')
%               After 30seconds simulation [#successfulLinks,UpperBound] = [100,138], after 5min [104,124.5]
%           params.Method = 4 ('Concurrent Deterministic')
%               After 30seconds simulation [#successfulLinks,UpperBound] = [100,138], after 5min [104,124.5]
%           params.ImproveStartTime doesn't make any difference
%          (2017/02/24): 
%               After clearing all insignificant values from 'A', params.Method = 1 ('Dual-simplex') seems to be the winning one. Time for achieving [106,138] reduced from 101s to 30s
%       3. (2017/02/14): [N,F,T] = [20,20,1]
%           Upon focusing more upon bound, the #successfulLinks reduces from 106 to 96 !!
%           Improving Parameters: method = 2
%           Unchanging Parameters: ImproveStartTime, Presolve = 2, SubMIPNodes = 1e6, MIPFocus=1, Heuristics = 1
%           Worsening Parameters: Cuts = 0; Symmetry
%       3. (2017/03/08): [N,F,T] = [20,20,1]
%           halfDuplexLTEacir: nSuccessfulLinks = 38 with 17% optimalty gap in 5min
%
%%%% -TODO:
%       1. Generalise function 'AXBproblem_to_Dxb_translate' to include F also (for the case F~=N)
%       2.  % -Scaling Up Noise variance, so that the (power values required) > (gurobi tolerance)  %% -TODO_Future  (Moved this to function 'optimalAlgorithmGurobi')
%           configDefault.noiseVariance = 1e5 * configDefault.noiseVariance;  %% - TODO_Urgent: Need to reduce this scaling
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = sch_gurobi(config,topology,txVUEPerRBPerTime_0,powerPerTxPerPerTime)


try

global probe;

N = config.N;
Ntx = config.N;
Nrx = config.N;
F = config.F;
T = config.T;
if(~gf.isempty('powerPerTxPerPerTime'))
    assert(size(powerPerTxPerPerTime,3)<=1,'Error txPower on multiple timeslots not supported at this moment');
    P = diag(powerPerTxPerPerTime);
else
    P = config.Pmax; 
end
PH = P*topology.H;

% % error('Sample error to test exception');

%% Schedule along time axis: Create list of VUEs scheduled on each time-slot

% Create problem Model
persistent Nold Fold Told model;
if(Nold==N & Fold==F & Told==T & ~gf.isempty('model'))   %% Made it this way, so that it can handle empty Nold, ... etc
    %% DO NOTHING
else    
    clear model;
    I = eye(F);
    A = (I-config.gammaT*topology.Lambda);  
    B = PH; 
    gamma = config.gammaT*config.noiseVariance;
    % Scale A,B, gamma
    Btemp = B; Btemp(logical(eye(size(B)))) = 0;
    scaleFactorA = 10^floor(log10(mean(abs(A(:)))));
    scaleFactorB = 10^floor(log10(mean(abs(Btemp(:)))));
    A = A/scaleFactorA;
    B = B/scaleFactorB;
    gamma = gamma/(scaleFactorA*scaleFactorB);
    % Create eta
%     eta = ceil(max(nansum(abs(A),2))*max(nansum(abs(B),1)));
    eta = gf.max(abs(A)*ones(size(A,2),size(B,1))*abs(B));
    eta = 10^ceil(log10(eta));      %% ceiling to the closest 10^x
    % -Create symbolic variables
    X = sym('X',[F,N,T]);
    Y1 = sym('Y1',[F,N,T]);
    % Split into 2 cases, when T=1 and T>1
    if(T==1) 
        % Warning: One catch here is that, if self-channel made too high value, then this algorithm also include self link as successful links, resulting in maximizing #VUEs to schedule.
        %       However not a severe problem whenver "F >= N", since all VUEs get scheduled normally
        A1 = -kron(B',A);  b1 = -gamma*ones(size(A1,1),1); 
        A1(logical(eye(size(A1)))) = 0;                                     %% This is to ignore self links while maximizing number of successful links, especially when self channel is nonzero.
        %if(~config.isEnableFullDuplex)                                      %% Ignore jth VUE reception link, if jth VUE is already scheduled in any of the RBs.
        %    A1 = A1 + eta*kron(eye(N)',ones(F));
        %end
        toleranceFactor = 1e-5; b1 = b1 - toleranceFactor*abs(b1);
        A2 = [A1; kron(eye(Ntx),ones(1,F)); kron(ones(1,Ntx),eye(F))];      %% Adding constraint for a VUE to be scheduled to max one time in a timeslot, and an RB can accommodate max 1 VUE
        b2 = [b1;ones(Ntx+F,1)];
        A3 = [A2,-eta*eye(size(A2,1),size(A1,1))];                          %% NOTE: To have the legacy implementation of the half duplex, multiply by 2*eta
        
        A = A3;                                                             %% i^th row of A1 correspnds to i%F RB and ceil(i,N)^th RX VUE
        b = b2;                                                             %% If j <= F*Nrx, then j^th column of A corresponds to j%F th RB and ceil(j,N) Rx VUE
        
        % Enforcing Half-Duplex
        if(~config.isEnableFullDuplex)                                      %% Ignore jth VUE reception link, if jth VUE is already scheduled in any of the RBs.
            A4 = zeros(N*F^2,size(A,2));                                    %% This is done by adding constraint X_{f,j} + Y_{fp,j} <= 1    for all {f,fp,j}.   (ie. eg 34 in the paper)
            b4 = zeros(N*F^2,1);                                            %%  Observe that 'Y' matrix in this code is complement of that of the paper
            r = 1;                                                        
            for j=1:Nrx
                for f=1:F
                    for fp=1:F
                        A4(r,F*(j-1)+f) = 1;
                        A4(r,F*Nrx + F*(j-1)+fp) = -1;                        
                        r = r+1;
                    end
                end
            end
            A = [A;A4];
            b = [b;b4];
        end

        c = zeros(size(A,2),1);
        c(F*Nrx+1:end) = 1;
        % Dummy creation of symbolic variables
        variablesALL = [X(:);Y1(:)];

        % -Remove insignificant values in 'A' (2017/02/24: This way, we are able increase sparsity of 'A' from 56% to 75%)
        toleranceOfInsignificantValuesSum = toleranceFactor*min(abs(b(b~=0)));
        [Asorted,AsortIndices] = sort(abs(A),2);
        isIgnore_Asorted = cumsum(Asorted,2) <= toleranceOfInsignificantValuesSum;
        isIgnore_A = nan(size(A));
        isIgnore_A(sub2ind(size(A),repmat([1:size(A,1)]',1,size(A,2)),AsortIndices)) = isIgnore_Asorted;
        A(logical(isIgnore_A)) = 0;
        % -Upperbound high values in A (Highly useful for the case config.isEnableFullDuplex = false )
        A = upperboundHighestValuesInA(A,b);
        
    else 
        assert(config.isEnableFullDuplex,'Error: Currently there is no provision to block reception (not necessarily self reception) of a scheduled VUE!!');
        % -Create symbolic variables
        Y2 = sym('Y2',[F,N,T]);
        V = sym('V',[F,N,N,T]);
        W = sym('Y1',[N,N]);
        variablesALL = [X(:);Y1(:);Y2(:);V(:);W(:)];        
        inequalities = [];  equalities = [];
        % 1a constraints  (ie, max number of RBs scheduled to a VUE, and max number of VUEs scheduled to an RB)
        inequalities = cat(1,inequalities, gf.vec(gf.sum(X,1)-1));
        inequalities = cat(1,inequalities, gf.vec(gf.sum(X,2)-1));
        % 1b constraints
        for t=1:T
            if(config.isEnableFullDuplex)
                inequalities = cat(1,inequalities, gf.vec(-A*X(:,:,t)*B - 2*eta*Y1(:,:,t) + eta*X(:,:,t) + gamma));                          %% (2017/03/02):  "eta*X(:,:,t)" is added to ignore self-links, especially when self channel is nonzero.
            else
                inequalities = cat(1,inequalities, gf.vec(-A*X(:,:,t)*B - 3*eta*Y1(:,:,t) + eta*X(:,:,t) + eta*ones(F,F)*X(:,:,t) + gamma)); %%                "eta*ones(F,F)*X(:,:,t)" is added to make it half duplex
            end
        end                                                                                                                                 
        % 2nd constraints
        equalities = cat(1,equalities,gf.vec(Y2-Y1));
        % 3rd constraints
        for t=1:T
            for j=1:N
                inequalities = cat(1,inequalities, gf.vec(-V(:,:,j,t) + X(:,:,t) + repmat(Y2(:,j,t),[1,N]) -1));
                inequalities = cat(1,inequalities, gf.vec(V(:,:,j,t) - X(:,:,t)));
                inequalities = cat(1,inequalities, gf.vec(V(:,:,j,t) - repmat(Y2(:,j,t),[1,N])));
            end
        end
        % 4th constraints
        inequalities = cat(1,inequalities, gf.vec(V - repmat(shiftdim(W,-1),[F,1,1,T])));
        inequalities = cat(1,inequalities, gf.vec(W - shiftdim(gf.sum(V,[1,4]),1)));
        c = -ismember(variablesALL,W(~eye(size(W))));        
        % Create A,b,c matrices
        A = gf.findCoefficientsPerVariable([inequalities;equalities],[variablesALL;1]);
        [A,b] = deal(A(:,1:end-1),-A(:,end));
        model.sense = [repmat('<',[numel(inequalities),1]); repmat('=',[numel(equalities),1]);];
    end
    model.A = sparse(A);
    model.rhs = b;
    model.obj = c;
    model.lb = zeros(size(A,2),1);
    model.ub = ones(size(A,2),1);
    model.vtype(1:size(A,2)) = 'B';
    Nold = N; Fold = F; Told = T;
end

if(~gf.isempty('txVUEPerRBPerTime_0'))
   model.start = generateStartingPoint(txVUEPerRBPerTime_0,variablesALL,X); 
end

% Solve the Problem
params = struct();
if(config.isEnableFullDuplex == false)
    params.FeasibilityTol = 1e-9;
    params.IntFeasTol = 1e-9;
else
    params.FeasibilityTol = 1e-6;
    params.IntFeasTol = 1e-6;    
end

% params.Method = 2;
params.SubMIPNodes = 1e6;   %% Stress more upon good quality feasible solution, rather than finding out bounds
params.MIPFocus=1;          %% Stress more upon good quality feasible solution
params.Heuristics = 1;      %% Stress more upon reducing objective value
params.TimeLimit = 6000;     %% disp(' change params.TimeLimit');
params.OutputFlag=0;        %% To turn off display
[x,out.exitFlag,result] = gurobi_(model,params);
% hasTimerRunOut = evalWithTimer('[x,out.exitFlag,result] = gurobi_(model,params);',params.TimeLimit);
Xv = reshape(x(ismember(variablesALL,X(:))),[F,N,T]);

Xv = gf.round(Xv,params.IntFeasTol);
assert(gf.all(sum(Xv,1)<=1),'Error: VUEs are scheduled to multiple RBs in a timeslot!!');
assert(gf.all(sum(Xv,2)<=1),'Error: Multiple VUEs scheduled to single RB!!');

% Assign the output
out.txVUEPerRBPerTime = zeros(F,T);
% if(hasTimerRunOut==0 && (out.exitFlag == 1 || out.exitFlag == 0))
if(out.exitFlag == 1 || out.exitFlag == 0)              %% NOTE: Here we are taking suboptimal solutions too...
    for t=1:T
        [rbIndices,txVUEs] = find(Xv(:,:,t));
        out.txVUEPerRBPerTime(rbIndices,t) = txVUEs;
    end
else
    out.txVUEPerRBPerTime = zeros(F,T);
end
    

catch ME
    ME = gf.saveWorkspaceUponException(ME);
    throw(ME);
end


end


%% Generate a starting point for Gurobi solver
%   TODO_Future: include starting vector for Y1,Y2,V & W too
function start = generateStartingPoint(txVUEPerRBPerTime_0,variablesALL,X)
    [F,N,T] = size(X);
    Xv = zeros(size(X));
    for t=1:T
        Xv(gf.sub2ind([F,N,T],[1:F]',txVUEPerRBPerTime_0(:,t),t)) = 1;
    end
    start = zeros(numel(variablesALL),1);
    start(ismember(variablesALL,X(:))) = Xv(:);
end



%% -Manipulate diagonal elements, so that diagonal element is not too big compared to other values in the row.
%    This to reduce the ill-condiitoning in the problem model in sch_gurobi(), when config.isEnableFullDuplex = false;
%   Warning: Here we assume "equal Power"
function H = upperBoundDiagonalValues(H)
    Htemp = H;
    Htemp(logical(eye(size(H)))) = 0;
    nonDiagonalSum = sum(abs(Htemp),2);
    nonDiagonalSum = 10.^ceil(log10(nonDiagonalSum));
    for iRow=1:size(H,1)
        if(abs(H(iRow,iRow))>nonDiagonalSum(iRow))
            H(iRow,iRow) = sign(H(iRow,iRow))*nonDiagonalSum(iRow);
        end
    end
end

%% -Pick each row, upperbound the highest value if it's above the dynamic range of other values in that row, and 'b'
%
function A = upperboundHighestValuesInA(A,b)
    B = abs(A);
    [maxValuePerRow,colIndexOfMaxValues] = max(B,[],2);
    for iRow = 1:size(A,1)
        row = B(iRow,:);
        row(colIndexOfMaxValues(iRow)) = 0;
        upperbound = sum(row)+abs(b(iRow));
        upperbound = 10^ceil(log10(upperbound));
        if(maxValuePerRow(iRow)>upperbound)
            A(iRow,colIndexOfMaxValues(iRow)) = sign(A(iRow,colIndexOfMaxValues(iRow)))*upperbound;
        end
    end
end



