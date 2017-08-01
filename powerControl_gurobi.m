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
%%%% -DEFINITIONS:
%       1. A link is a physical link present for all RBs to all receivers in all timeslots. Hence, there are 'F x Nrx x T' currentLinks in total.
%           lth link means the transmissionfrom on 'l%F' th RB to 'floor(l,Nrx)th' VUE. Note that this doesn't tell anything about transmitter VUE
%
%%%% -OPTIONS:
%
%%%% -EXAMPLES:
%       1. 
%
%%%% -NOTES:
%
%%%% -NOTES (Programming):
%
%
%%% -Exection
%       1.
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


function powerPerTxPerPerTime = powerControl_gurobi(config,topology,txVUEPerRBPerTime,powerPerTxPerPerTime)

toleranceFactor = 1e-5;     %% To accommodate minor roudoff errors in Gurobi, we are making constraints more stricter apriori :-)

Ntx = config.N;
Nrx = config.N;
F = config.F;
T = config.T;

if(~gf.any(txVUEPerRBPerTime))
    powerPerTxPerPerTime = nan(Ntx,1,T);
    return;
end

[A,b] = topology.computeConstraintMatrixForPowerValues(txVUEPerRBPerTime);
assert(size(A,3)==T,'Error: 3rd dimension of ''A'' must be along timeslots!!');

% -Make A to 2D matrix (ie. if T=2, then make A = [A(:,:,1),0; 0,A(:,:,2)])
indexMatrix = logical(kron(eye(T),ones(size(A(:,:,1)))));
Anew = zeros(size(indexMatrix));
Anew(indexMatrix) = A(:);
A = Anew;
b = b(:);

% -Get all parameters related to currentLinks
lS = topology.getLinkStatistics(txVUEPerRBPerTime,A);
%A = A(lS.linksToConsider,:);
%b = b(lS.linksToConsider);

% -Scale 'A' & 'b' so as to get rid off rouding issues in Gurubi
scalingFactor = (1e-3)/mean(abs(b));
scalingFactor = 10^ceil(log10(scalingFactor));      %% ceiling to the closest 10^x
A = scalingFactor*A;
b = scalingFactor*b;

% -Prepare for gurobi
params = struct();
params.FeasibilityTol = 1e-9;      %% Necessary since dynamic range b/w A & b is 10^10 !!
params.IntFeasTol = 1e-9;
params.TimeLimit = 600;
params.OutputFlag=0;      %% To turn off display place '0', otherwise '1'

% -Compute 'powerPerTxPerPerTime'  (NOTE 2017/04/17: Optimising on each timeslot, actually gives better results)
powerPerTxPerPerTime = zeros(Ntx,1,T);
nLinks = numel(lS.txVUEPerLink);
isSuccessPerLink = false(nLinks,1);
for t=1:T
    currentLinks = intersect(lS.linksToConsider,find(lS.timeslotPerLink==t));
    if(~isempty(currentLinks))
        A_ = A(currentLinks,(t-1)*Ntx+1:t*Ntx);
        beta = 1 / (100*Ntx*T*config.Pmax);
        beta = 10^floor(log10(beta));      %% Flooring to the closest 10^x
        c = beta*ones(1,size(A_,2));
%         c = [];
        [x,~,exitflag,result] = gf.linprog(c,A_,b(currentLinks)-toleranceFactor*abs(b(currentLinks)),[],[],zeros(Ntx,1),config.Pmax*ones(Ntx,1),'-gurobi','-params',params);
        if(exitflag==1)
            powerPerTxPerPerTime(:,1,t) = x;
            % -Find out successful links, and remove those from 'lS.linksToConsider'
            successfulTx = lS.txVUEPerLink(currentLinks(result.x(Ntx+1:end)==0));
            successfulRx = lS.rxVUEPerLink(currentLinks(result.x(Ntx+1:end)==0));
            for i=1:numel(successfulTx)
                isSuccessPerLink(lS.txVUEPerLink==successfulTx(i) & lS.rxVUEPerLink==successfulRx(i)) = 1;
            end        
            lS.linksToConsider = setdiff(lS.linksToConsider,topology.convertRealLinksToLinks(lS,lS.realLinkPerLink(isSuccessPerLink)));
        end
    end
end

%% -Debug:
% isSuccessPerTxPerRx = topology.findSuccessStatus(lS,isSuccessPerLink);
% disp(['powerControl_gurobi: ',num2str(gf.nansum(isSuccessPerTxPerRx))]);
% isSuccessPerTxPerRx


end