%%%%%%%%%%%%%%%%%% -COPYRIGHTS & WARRANTY- %%%%%%%%%%%%%%%%%%%%%%%%
%% It is provided without any warranty of fitness
%% for any purpose. You can redistribute this file
%% and/or modify it under the terms of the GNU
%% Lesser General Public License (LGPL) as published
%% by the Free Software Foundation, either version 3
%% of the License or (at your option) any later version.
%% (see http://www.opensource.org/licenses for more info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% -NOTES
%   1. txVUE/rxVUE = 0   => No one is transmitting/receiving 
%   2. 'config' for system model, and 'opt' for user inputted options (or parameters)
%
% -Execution
%   1. opt.Fs = [5,10,20]';  opt.Ts = [1,2,3,4]';  opt.Ns = 20; opt.nTrials = 50;    =>  2.5 hours
%   2. opt.Fs = [5,10,20]';  opt.Ts = [1,2,3,4,5]';  opt.Ns = 20; opt.nTrials = 50;  =>  3.8 hours
%   3. 'sch_gurobi' is not finishing even after one hour for [N,F,T] = [20,20,1], nTrials = 1 !!!!
%   4. opt = struct('nTrials',100,'Ns',[3:20]','Fs',[5,10,20]','Ts',[1:5]', 'schAlgorithmsToRun', {{'sch_naturalScheduling','sch_blockInterleaver1','sch_blockInterleaver2','sch_blockInterleaver3','sch_Heuristic1'}}) 
%           => 124 hours
%       (Update 2017/02/09): => 42 hours for a single core machine
%   5. out=architecture2_wrapper([],struct('nTrials',1,'Ns',[20]','Fs',[5,10,20]','Ts',[1:5]', 'schAlgorithmsToRun', {{'sch_naturalScheduling','sch_Heuristic1','sch_Heuristic1([],[],[],[],true)'}},  'powerControlAlgorithmsToRun',{{'powerControl_kangWang'}}));
%       (2017/02/28): => 166 seconds
%   6. out1=architecture2_wrapper('',struct('nTrials',100))
%      gf.nansum(nanmean(gf.congregate('out1{:}(:).isSuccessPerTxPerRxPerAlgorithm'),[2]),[3,4])
%
% -TODO
%
%% -Author
%   Anver Hisham <anverhisham@gmail.com>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = architecture2_wrapper(config,opt)

tic;
addpath('~/Programming/Matlab/MyMatlabFunctions/');
addpath('~/Programming/Matlab/ALLtoALL_Communication/Library');
warning('off','MATLAB:singularMatrix');         %% Disable printing warnings
warning('off','MATLAB:nearlySingularMatrix');  

clearvars -global probe
global probe;
probe.isDebug = false;

% Glenn inputs
isGlenn = isunix && ~ismac;
% isGlenn = true; disp('#DeleteGlennTesting isGlenn made to true');       %% #DeleteGlennTesting
if(isGlenn)
   if(gf.isempty('opt.Ns')) opt.Ns = [3:20]'; end
   if(gf.isempty('opt.Fs')) opt.Fs = [6,12,20]'; end
   if(strcmpi(opt.schAlgorithmsToRun{1},'sch_gurobi'))
    if(gf.isempty('opt.Ts')) opt.Ts = [1]'; end
    if(gf.isempty('opt.nTrials')) opt.nTrials = 1; end
   else
    if(gf.isempty('opt.Ts')) opt.Ts = [1:20]'; end              %% (2017/7/13) #Log: Changed from [1:5]' to [1.20]'
    if(gf.isempty('opt.nTrials')) opt.nTrials = 100; end
   end
end
opt = patchOpt('opt',isGlenn);
% opt.nThreads=1; opt.nTrials = 1; disp('#DeleteGlennTesting nTrials made to 1');         %% #DeleteGlennTesting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

config = patchConfig('config');
if(opt.nThreads>1) probe.isDebug = false; end  

in = createInputForexecuteIteratively(config,opt);
in = nullifyNonMarginalConfigurations(in);      %% We are interested in results of only marginal configurations, ie F=20,T=1,N=20

%% -Simulate. (Resulting out is "<#VUEs> <#Trials> <#Algorithms>")
% try
% persistent architecture2_wrapper_tryCatchDisplay; 
% if(isempty(architecture2_wrapper_tryCatchDisplay)) architecture2_wrapper_tryCatchDisplay = true;  disp('Warning: No try catch block in architecture2_wrapper!!'); end

out = gf.executeIteratively(@architecture2,in,opt.nThreads,'-reverse');
[uniqueID,summary] = getUniqueIDandSummary(config,opt,'architecture2_wrapper.log');
disp(summary);
save([uniqueID,'_architecture2_wrapper.mat'],'out','config','opt','probe','-v7.3');

% catch ME
%    ME.getReport 
%    toc
%    rethrow(ME);
% end

end


function in = nullifyNonMarginalConfigurations(in)
    for f=1:size(in,1)
        for t=1:size(in,2)
            for n=1:size(in,3)
                if(double(f==size(in,1)) + double(t==1) + double(n==size(in,3)) < 2)    %% If atleast 2 conditions are satisfied, then don't reset.
                    in{f,t,n}{1}.F = 0;
                    in{f,t,n}{1}.T = 0;
                    in{f,t,n}{1}.N = 0;
                end
            end
        end 
    end
end



%% -Make 'opt' by varying #RBs, #Time-slots and #VUEs 
function in = createInputForexecuteIteratively(config,opt)
    in = cell(numel(opt.Fs),numel(opt.Ts),numel(opt.Ns));
    optCallee = rmfield(opt,{'Ns','Fs','Ts','nThreads'});
    [in{:}] = deal({config,optCallee});
    for f=1:numel(opt.Fs)
        for t=1:numel(opt.Ts)
            for n=1:numel(opt.Ns)
                in{f,t,n}{1}.F = opt.Fs(f);
                in{f,t,n}{1}.T = opt.Ts(t);
                in{f,t,n}{1}.N = opt.Ns(n);
            end    
        end 
    end
end

%% -
function [uniqueID,summary] = getUniqueIDandSummary(config,opt,fileName)
    uniqueID = datestr(now,'yyyymmddHHMMssFFF');
    summary = ['   ',10,10,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',10];
    summary = [summary,'uniqueID = ',uniqueID,10];
    summary = [summary,gf.struct2string(opt,true)];
    summary = [summary,gf.struct2string(config,true)];
    if(~gf.isempty('config.configChannel'))
        summary = [summary,gf.struct2string(config.configChannel,true)];
    end
    summary = [summary,evalc('toc')];   
    if(~gf.isempty('fileName'))
        fid = fopen(fileName,'a');
        fprintf(fid,summary);
        fclose(fid);
    end
end


function opt = patchOpt(optString,isGlenn)
    if(evalin('caller',['gf.isempty(''',optString,''')']))
        opt = struct();
    else
        opt = evalin('caller',optString);
    end;    
    optDefault.nTrials = 1; 

    optDefault.schAlgorithmsToRun = {'sch_naturalScheduling'};                             %% -MAUAL INPUTS %%%%%%%%%%%%%%%%%%%%
    optDefault.powerControlAlgorithmsToRun =  {'powerControl_gurobi'};
    optDefault.Ns = [3,4,5,6,7]';  optDefault.Fs = [20]';  optDefault.Ts = [1]';   %% Look for #AffectedBy_opt.Fs, #AffectedBy_opt.Ts while changing optDefault.Fs
    
    if(isGlenn) optDefault.nThreads=16; else optDefault.nThreads=1; end             %% For Glenn cluster
    opt = gf.mergeTwoStructures(optDefault,opt,'',true);
end

%% -
function config = patchConfig(configString)
    if(evalin('caller',['gf.isempty(''',configString,''')']))
        config = struct();
    else
        config = evalin('caller',configString);
    end;
    %% -Default Parameters
    configDefault.N = 4;
    configDefault.F = 2;
    configDefault.T = 3;
    configDefault.gammaT = db2pow(5);
    configDefault.acirType = 'LTEmask';         %'scFDMA';
    configDefault.isEnableFullDuplex = false;        %%
    %% -Merging Inputted config with Default inputs,
    config = gf.mergeTwoStructures(configDefault,config,'',true);
end
