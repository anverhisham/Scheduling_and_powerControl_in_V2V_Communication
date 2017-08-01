%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% -Inputs 
N = 5;      %%  #VUEs
F = 5;      %%  #Resource Blocks
T = 1;      %%  #Subframes (or Total time period in milli seconds) 
FT = F*T;
gammaT = sym('gammaT');
I = eye(N);

P = diag(sym('P',[N,1]));
H = sym('H%d%d',[N,N]);
alphaB = sym('a%d%d',[F,F]);
% alphaB(logical(eye(size(alphaB)))) = 0;
str = repmat(['alphaB,'],[1,T]);
str(end) = [];

alpha = eval(['blkdiag(',str,');']);
% X = sym('X',[FT,N]);
X = eye([FT,N]); X = X + X([end,1:end-1],:);    %% Scheduling is: VUE5, VUE1, 2, 3, 4
% X = [1,0,0;1,0,0; 0,1,0;0,1,0;  0,0,1;0,0,1;];

% LHS = P*H - gammaT*X.'*alpha.'*X* P*H;
% LHS = P*H - X.'*alpha.'*X* P*H;
LHS = X*P*H - alpha.'*X* P*H;

LHS(1,3)


