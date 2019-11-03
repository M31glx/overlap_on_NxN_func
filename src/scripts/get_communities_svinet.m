% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% setup call

binpath = '/usr/local/bin/svinet' ;
    
%% the svinet opts  
% SVINET: fast stochastic variational inference of undirected networks
% svinet [OPTIONS]
% 	-help		usage
% 
% 	-file <name>	input tab-separated file with a list of undirected links
% 	-n <N>		number of nodes in network
% 	-k <K>		number of communities
% 	-batch		run batch variational inference
% 	-stratified	use stratified sampling
% 	 * use with rpair or rnode options
% 	-rnode		inference using random node sampling
% 	-rpair		inference using random pair sampling
% 	-load-validation <fname>		use the pairs in the file as the validation set for convergence
% 	-load-test <fname>		use the pairs in the file as the test set
% 	-label		tag output directory
% 	-link-sampling	inference using link sampling 
% 	-infset		inference using informative set sampling
% 	-preprocess	preprocess to run informative set sampling
% 	-rfreq		set the frequency at which
% 	 * convergence is estimated
% 	 * statistics, e.g., heldout likelihood are computed
% 	-max-iterations		maximum number of iterations (use with -no-stop to avoid stopping in an earlier iteration)
% 	-no-stop		disable stopping criteria
% 	-seed		set GSL random generator seed
% 	-gml		generate a GML format file that visualizes link communities

maxiters = 500 ;
ncommunities = 10 ;

%% run svinet and recover communities

rng(123)

dat_thr = 0.10 ;

% the file we want
inputfile = [ OUTDIR_PROC '/groupavg_' PARC_STR ... 
        '_edgelist_thr' sprintf('%0.2f',dat_thr) '.txt' ] ;

exe_call = strcat(binpath, ...
    ' -file', 32, inputfile, ...
    ' -n', 32, num2str(NUM_NODES), ...
    ' -k', 32, num2str(ncommunities), ...
    ' -max-iterations', 32, num2str(maxiters),...
    ' -link-sampling -no-stop ' ) ;
    
numRun = 100 ;
svinet_groups = cell(numRun,1) ;
svinet_comms = cell(numRun,1) ;

for idx = 1:numRun 
    disp(idx)
    
    tmpDir = strcat(OUTDIR_INTERM,'/tempDir/') ;
    mkdir(tmpDir) 
    cd(tmpDir)

    % run the command
    [~,~] = system([ exe_call ' -seed ' num2str(randi(1000)) ]) ;

    outDir = dir([tmpDir '/n*']) ; 
    %outDir(1:2) = [] ;
    
    % groups data
    outGroupPath = [ outDir.folder '/' outDir.name '/groups.txt' ] ;
    outGroupTable = readtable(outGroupPath) ;
    outGroup = table2array(outGroupTable(:,2:end)) ;
    [~,sss] = sort(outGroup(:,1)) ;
    outGroup = outGroup(sss,2:end) ;
    svinet_groups{idx} = outGroup ;

    % communities data
%     outCommPath = [ outDir.folder '/' outDir.name '/communities.txt' ] ;
% %     outComm = table2cell(readtable(outCommPath,'ReadVariableNames',0)) ;
% %     outComm = cellfun(@(a)str2num(a),outComm,'UniformOutput',false) ;
% %     outComm = comms_cell_2_mat(outComm,NUM_NODES) ;
%     outComm = dlmread(outCommPath) ;
%     outComm = comms_mat_2_mat(outComm,NUM_NODES) ;
%     svinet_comms{idx} = outComm ;     
    svinet_comms{idx} = outGroup > COMM_PROB_THR ;    

    cd(PROJECT_DIR)    
    rmdir(tmpDir,'s')

end

%% now get the centroid using 

% make a temp dir
tempdir = strcat(OUTDIR_INTERM,'/tempcomms/') ;
mkdir(tempdir)

% first write out all the files
for idx = 1:numRun
    disp(idx)
    write_comms(svinet_comms{idx},[ tempdir '/comm' num2str(idx) '.txt'])
end

% now compare all communities
ovrmutinfomat = zeros(numRun) ;
mutexe = [ PROJECT_DIR '/src/external/mutual3/mutual' ] ;
for idx = 1:numRun
   for jdx = 1:numRun
       if idx >= jdx
           continue
       else
           disp([num2str(idx) ' ' num2str(jdx)])
       end
       
       file1 = [ tempdir '/comm' num2str(idx) '.txt'] ;
       file2 = [ tempdir '/comm' num2str(jdx) '.txt'] ;      
       
       [a,b] = system([ mutexe ' ' file1 ' ' file2 ' | awk  ''{print $2}'' ']) ;
       ovrmutinfomat(idx,jdx) = str2double(b) ;
   end
end

rmdir(tempdir,'s')

%  find the centroid
simmat = ovrmutinfomat + ovrmutinfomat' ;
[~,centind] = max(sum(simmat)) ;

svinet_cent = svinet_comms{centind} ;

%% try to algin all the comms

% didn't work, doesn't look good...

% [~,svinet_cent_greedy] = max(svinet_groups{maxind},[],2) ;
% 
% svinet_groups_align = cell(numRun,1) ;
% % align all comms to the cent_greedy
% for idx = 1:numRun
%     disp(idx)
%     
%     [~,tmpGreedyComms] = max(svinet_groups{idx},[],2) ;
%     [~,reorder] = hungarianMatch(svinet_cent_greedy,tmpGreedyComms) ;
%     
%     svinet_groups_align{idx} = svinet_groups{idx}(:,reorder) ;
% end

%% look at it
 
% for idx = 1:numRun
% imagesc(svinet_groups_align{idx})
% waitforbuttonpress
% end

%% make an agreement

agree_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;
agreew_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;

for idx = 1:numRun
    disp(idx)
    
    agreew_dat(:,:,idx) = svinet_comms{idx} * svinet_comms{idx}' ;
    agree_dat(:,:,idx) = int32(agreew_dat(:,:,idx) > 0) ;

end

%% save results

save([ OUTDIR_PROC '/svinet_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ],...
    'svinet_comms','svinet_groups','svinet_cent','centind') ;



