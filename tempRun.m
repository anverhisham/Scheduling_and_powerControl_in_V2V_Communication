%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% foldername = './Results/2017March03/fullDuplexLTEacir/interleaver1/';
% analyzeResult.sch_interleavers_widths(foldername);
% 
% foldername = './Results/2017March03/fullDuplexLTEacir/interleaver2/';
% analyzeResult.sch_interleavers_widths(foldername);
% 
% foldername = './Results/2017March03/fullDuplexLTEacir/interleaver3/';
% analyzeResult.sch_interleavers_widths(foldername);

% % tic
% % architecture2_wrapper(struct('acirType','scFDMA','isEnableFullDuplex',false),struct('nTrials',1,'Ns',[20]','Fs',[20]','Ts',[1]', ...
% %     'schAlgorithmsToRun', {{'sch_naturalScheduling','sch_Heuristic1','sch_Heuristic1([],[],[],[],2)'}},  'powerControlAlgorithmsToRun',{{'powerControl_kangWang'}}))
% % toc



tic, analyzeResult.majorAlgorithms('~/Programming/Matlab/ALLtoALL_Communication/Architecture2/Results/2017March03/halfDuplexLTEacir/'); toc

tic, out_halfDuplexLTEacir=architecture2('',struct('nTrials',10)); toc, save('out_halfDuplexLTEacir.mat','out_halfDuplexLTEacir');