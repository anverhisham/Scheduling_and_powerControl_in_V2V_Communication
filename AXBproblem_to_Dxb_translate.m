%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -This function translate the following problem
%       "compute best permutation matrix 'closestPermutationMatrix' s.t. |AXB >= gamma| is maximized"
%   into
%       "min c'x  s.t.  Dx<=b_C" where x(1:end/2) is permutations matrix X, and x(end/2:end) is infeasiblity vector y (ie yi=1 => i^th link is infeasible ) 
%       Note that extra bolean constraint for 'x' has to be placed by the client program (what does this mean?)
%       Note that the returned constraint set [D,b] are normalised w.r.t mean absolute value 
%
% translationType = 0       c = [0,0,.., 1,1,...]
%  D = | -kron(B',A),            eta2*I(N^2,N^2)  |     b = | -gamma*1(N^2,1)  | 
%      | col & row sum inequality,   0(2N,N^2)   |         |        1(2N,1)   | 
%      | -total sum inequality,      0(1,N^2)    |         |            -N         | 
%      | -Positivity inequality,     0(N^2,N^2)  |         |        0(N^2,1)       | 
%
% translationType = 1
%  D = | -kron(B',A),              |     b =| -gamma*1(N^2,1) | 
%      | col & row sum inequality |         |        1(2N,1)  | 
%      | -total sum inequality    |         |            -N   | 
%      | -Positivity inequality   |         |        0(N^2,1) | 
%
% translationType = 2
%       AXB >= X(gamma-eta*(1-Yc))
%   =>  (kron(B',A)+kron(Zceta,I)) X(:) >= gamma     {Zceta = eta*(1-Yc)}  
%       Then remove eta containing constraints, and add Birkhoff polytope. 
%
% -NOTE
%   1. Self links (VUE i to VUE i) are always removed in this function.
%
% -Warning:
%   The structure of 'D' matrix shouldn't be disturbed. (ie Top-left quarter of 'D' must corresponds to link success constraint)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [c,D,b] = AXBproblem_to_Dxb_translate(A,B,gamma,Yc,translationType)

    N = size(A,1);
    I = eye(size(A));
    if(gf.isempty('Yc')) Yc = ones(N); Yc(logical(eye(N))) = 0; end
    if(gf.isempty('translationType')) translationType = 0; end
    if(numel(gamma)==1) gamma = gamma*ones(N); end
    assert(all(size(A)==[N,N]) && all(size(B)==[N,N]) && all(size(gamma)==[N,N]),'Error(AXB_Algorithm1): Inputs must be square matrices');
    Yc = reshape(Yc(1:N^2),N,N);
    X = nan(N,N);
    
    %% -Converting the problem |AXB >= gamma| into Cx <= b_C
    C1 = -kron(transpose(B),A);
    b_C = -gamma(:);
    %%%% -STATUS: Now "AXB >= gamma" is converted into "Cx <= b_C"  %%%
    %% -Integrate Yc into A & B 
    if(translationType==2) eta = nan; else eta = 1e3*ceil(max(nansum(abs(C1),2)+abs(b_C))); end
    Zc = ~Yc; Zceta = eta*Zc;
    C2 = C1-(kron(Zceta',I));
    
    %% -Create Birkhoff Polytope set 'L' i.e.  L{1}x <= L{2}, and make to Dx <= b (NOTE(2016/10/13): Made it relaxed Birkhof Polytope, by relaxing the summation constraint) 
    L{1} = [gf.getRowSumMatrix(size(X)); gf.getColSumMatrix(size(X)); -ones(1,numel(X))];   %% -last inequality is for avoiding the case of highly sparse 'X'
    L{2} = [ones(size(X,1),1); ones(size(X,2),1); -N];                                    %%TODO_Test: Relaxed total sum inequality from '-N' to '-0.9*N' 
    % Create Positivity constraints Pc{1} and Pc{2} 
    Pc{1} = -eye(N^2);
    Pc{2} = zeros(N^2,1);
    %% -Append both Birkhoff Polytope constraint and Positivity constraints
    D = [C2;L{1};Pc{1}];   b = [b_C;L{2};Pc{2}];
    [D,b] = sf.normaliseConstraints(D,b);
    
    %% -Create objective vector c (in order to make it as a minimization problem cT*x)
    c = [];
    
    if(translationType==0)                      %% Addind off "maximize #satisfied constraints" part (Adding eta2 to the terms)
        c = zeros(size(D,2),1);
        eta2 = ceil(max(nansum(abs(D(1:N^2,:)),2)));
        %% -Converting "Cx <= b_C" to "[C1,-eta*I][x,y] <= b_C". 'y' is infeasiblity factor (y=1 => that link is infeasible)
        D = gf.cat(2,D,-eta2*eye(size(A,1)*size(B,2)),'-appendValue',0);  %% TODO_Urgent     
        c = [c;ones(size(D,2)-numel(c),1)];
    elseif(translationType==1)                  
    elseif(translationType==2)
        noncandidateLinks = find(any(isnan(D),2));
        D(noncandidateLinks,:) = [];
        b(noncandidateLinks) = [];
    else
        error('Provide a valid translationType!!');
    end
end



