%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% -Got this code from http://www.ece.uvic.ca/~wslu/Talk/SeDuMi-Remarks.pdf#page=7 
%
%  Solve Problem
%
%      min c'y
%
%      s.t 
%       Ay <= b
%       F{1} + y(1)F{2} + y(2)F{3} ... y(n)F{n+1} >= 0  (ie positive semidefinite!)  
%
% -Example-1
%   F = { -[2 -0.5 -0.6; -0.5 2 0.4; -0.6 0.4 3],  -[0 1 0; 1 0 0; 0 0 0],   -[0 0 1; 0 0 0; 1 0 0],   -[0 0 0; 0 0 1; 0 1 0], eye(3) };  
%   c = [0,0,0,1]';
%   [y, info] = solveSDPstandard(F,c);
%
% -Example-2
%   A = -[1 0 0 0; -1 0 0 0; 0 1 0 0; 0 -1 0 0; 0 0 1 0];  b = -[0.7 -1 0 -0.3 0]'; 
%   F = { -[2 -0.5 -0.6; -0.5 2 0.4; -0.6 0.4 3],  -[0 1 0; 1 0 0; 0 0 0],   -[0 0 1; 0 0 0; 1 0 0],   -[0 0 0; 0 0 1; 0 1 0], eye(3) };  
%   c = [0,0,0,1]';
%   [y, info] = solveSDPstandard(F,c,A,b);
%
% -Example-3
%   m=10; F{1} = randn(m,m); F{2} = randn(m,m); F{3} = randn(m,m);   c = [0,-1];   
%   [y, info] = solveSDPstandard(F,c);
%
%
% -PREREQUISITE:
%       Install SeDuMi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y, info] = solveSDPstandard(F,c,A,b)
if(gf.isempty('A')) A = []; end
if(gf.isempty('b')) b = []; end
p = length(c);
bt = -c;
ct = [b;vec(F{1})];
At = nan(numel(F{1}),p);
for i = 1:p
    At(:,i) = -gf.vec(F{i+1}); 
end
At = [A;At];
K.l = size(A,1);
K.s = size(F{1},1);
[~, y, info] = sedumi(At,bt,ct,K,struct('fid',0)); 

end



