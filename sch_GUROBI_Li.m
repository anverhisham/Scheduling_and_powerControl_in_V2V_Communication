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
%       1. Need to investigate why following example not showing proper result ???
%           out = sch_GUROBI_Li(struct('T',2,'F',3,'N',6)) 
%
%%%% -Versions
% 
%
%%%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: varargin is placed so that this function won't throw error even while calling with more inputs 
function out = sch_GUROBI_Li(config,~,~,~,blockInterleaverWidth,varargin)

if(gf.isempty('blockInterleaverWidth'))
   blockInterleaverWidth = 3; 
end
    
persistent Nold Fold Told model;
N = config.N;
F = config.F;
T = config.T;

%% Schedule along time axis: Create list of VUEs scheduled on each time-slot

% Create problem Model
if(Nold==N & Fold==F & Told==T)   %% Made it this way, so that it can handle empty Nold, ... etc
    %% DO NOTHING
else
    w = sym('w',[N,T]);         %% w(i,t) = 1 if ith VUE is scheduled on t^th time-slot 
    y = sym('y',[N,N,T]);       %% y(i,j,t) = 1 if ith VUE can communicate with jth VUE, on t^th time-slot 
    Y = sym('Y',[N,N]);         %% Y(i,j) = 1 if ith VUE can communicate with jth VUE, on any time-slot 
    variablesALL = [w(:);y(:);Y(:)];
    % Create an array of expressions for inequalities inequalities(:)<=0
    % Constraints on w_it
    inequalities = sum(w,1).'-F;
    % Constraints on y_ijt
    for i=1:N
        for j=1:N
            inequalities = cat(1,inequalities, shiftdim(y(i,j,:),2) -w(i,:).'-w(j,:).');
            inequalities = cat(1,inequalities, shiftdim(y(i,j,:),2) +w(i,:).'+w(j,:).' -2);
        end
    end
    % Constraints on Y_ij
    for i=1:N
        for j=1:N
            inequalities = [inequalities; Y(i,j)-sum(shiftdim(y(i,j,:),2))];
            inequalities = [inequalities; -Y(i,j)+shiftdim(y(i,j,:),2)];
        end
    end
    % Create A,b,c matrices
    A = gf.findCoefficientsPerVariable(inequalities,[variablesALL;1]);
    [A,b] = deal(A(:,1:end-1),-A(:,end));
    c = -ismember(variablesALL,Y(~eye(size(Y))));
    clear model;
    model.A = sparse(A);
    model.rhs = b;
    model.obj = c;
    model.lb = zeros(size(A,2),1);
    model.ub = ones(size(A,2),1);
    model.vtype(ismember(variablesALL,w(:))) = 'B';     %% TODO: verify this
    model.vtype(~ismember(variablesALL,w(:))) = 'C';
    Nold = N; Fold = F; Told = T;
end

% Solve the Problem
[x,exitFlag] = gurobi_(model);
w = reshape(x(1:N*T),N,T);
for t=1:T
    txVUEsPerTime{t} = find(w(:,t));
end

%% Schedule along frequency axis: by Block interleaving
txVUEsPerTime = cellfun(@(x) sort(x),txVUEsPerTime,'UniformOutput',false);                                      %% Sort on every time-slot
txVUEsPerTime = cellfun(@(x) gf.interleaveBlock_old(x,blockInterleaverWidth),txVUEsPerTime,'UniformOutput',false);  %% Do block interleaving
out.txVUEPerRBPerTime = cell2mat(txVUEsPerTime);

end