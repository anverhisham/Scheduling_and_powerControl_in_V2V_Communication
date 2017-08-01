%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -Find the max-volume ellipsoid (Eu+c, ||u||<=1) inside a Polytope Ax<=b
%       Got this code from https://groups.google.com/forum/#!topic/mpt-user/pQtJi34dXu4
%
% -NOTES:
%   1. To find min-volume ellipsoid containtaining Polytope we need set of vertices, not inequalities alone.
%       whereas inscribing ellipsoid can be found out either from vertices or from inequalities! 
% -EXAMPLES
%   1. Test example,
%       n=4^2; A=[eye(n);-eye(n)]; b=[ones(n,1);zeros(n,1)]; [c,E] = findMaxVolumeEllipsoidInPolytope_yalmip(A,b)
%       n=4^2; A=[eye(n);-eye(n)]; b=[ones(n,1);zeros(n,1)]; H = -2*eye(n); f = zeros(n,1); [x, fval, exitflag, output] = cplexqp (H, f, A, b) 
%
% -OBSERVATION:
%   1. This code takes hours of simulation even for 100 dimensions!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [c,E] = findMaxVolumeEllipsoidInPolytope_yalmip(A,b)

% find the largest inner ellipsoidal approximation. Here, the ellipsoid is 
% given by { Ex+c | ||x||_2 \le 1 } 
% (see Boyd, Vandenberge: Convex Optimization, 2004, section 8.4.2) 
    dim = size(A, 2);
    E = sdpvar(dim, dim); 
    c = sdpvar(dim, 1);       %% center-point of ellipse

    con = []; 
    for i = 1:size(A, 1) 
        con = con + [ norm(E*A(i, :)') + A(i, :)*c <= b(i) ]; 
    end 
    con = con + [ E>=0 ];
    diagonostic = optimize(con, -logdet(E),sdpsettings('verbose',0)); 

% find the smallest outer ellipsoidal approximation, i.e., find E and c 
% % % such that the ellipsoid { x | ||Ex+c||_2 \le 1 } contains all vertices of 
% % % the polytope 
% % % (see Boyd, Vandenberge: Convex Optimization, 2004, section 8.4.1) 
% % dim = size(A, 2);
% % E = sdpvar(dim, dim); 
% % c = sdpvar(dim, 1); 
% % con = []; 
% % for i = 1:size(A, 1)
% %     con = con + [ norm(E*P.V(i, :)'+c) <= 1 ]; 
% % end 
% % con = con + [ E>=0 ]; 
% % solvesdp(con, -logdet(E)); 

    E = double(E);
    c = double(c);

end


