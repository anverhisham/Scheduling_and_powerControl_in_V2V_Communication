%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tic;
addpath('~/Programming/Matlab/MyMatlabFunctions/');
addpath('~/Programming/Matlab/ALLtoALL_Communication/Library');
warning('off','MATLAB:singularMatrix');         %% Disable printing warnings
warning('off','MATLAB:nearlySingularMatrix');  

folder_mat = './';
% file_mat = 'temp_computeNSuccessfullLinks.mat';
file_mat = '20170227173304881_architecture2_wrapper.mat';
temp = load([folder_mat,file_mat]);
out = temp.out; opt = temp.opt; config = temp.config;

nSuccessPerFsPerTsPerNsPerAlgorithm = permute(nanmean(gf.nansum(gf.congregate('out{:,:,:}(:).isSuccessPerTxPerRxPerAlgorithm')==1,[5,6]),4),[1,2,3,7,4,5,6]); 
nSuccessPerFsPerTsPerNsPerAlgorithm = bsxfun(@rdivide,nSuccessPerFsPerTsPerNsPerAlgorithm,shiftdim(opt.Ns,-2));
nSuccessPerFsPerTsPerNsPerAlgorithm_ = matrixGeneralized([],'','',nan);
nSuccessPerFsPerTsPerNsPerAlgorithm_(opt.Fs,opt.Ts,opt.Ns,:) = nSuccessPerFsPerTsPerNsPerAlgorithm;

uniqueID = regexp(file_mat,'^\d*','match');
save([uniqueID{:},'_nSuccessPerFsPerTsPerNsPerAlgorithm.mat'],'nSuccessPerFsPerTsPerNsPerAlgorithm');

disp(['Execution of ',mfilename,' over ...']);

toc
