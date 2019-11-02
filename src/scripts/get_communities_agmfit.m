% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% setup call

binpath = '/home/jfaskowi/JOSHSTUFF/software/snap/examples/agmfit/agmfitmain' ;
    
%% the agmfit opts  
% ========================================================================
%     Community detection by Community-Affiliation Graph Model 
% ========================================================================
% 
% The example implements overlapping community detection by the 
% Community-Affiliation Graph Model (AGM).
% 
% This program detects network communities from a given network by fitting 
% AGM, a probabilistic generative model for networks, to the given network 
% by maximum likelihood estimation. User can specify how many communities 
% she would detect, or let the program automatically determine the number 
% of communities in the network based on the structure of the network.
% 
% Fitting procedure and the community-Affiliation Graph Model are 
% described in the following paper:
% J. Yang and J. Leskovec, Structure and Overlaps of Communities in 
% Networks, SNAKDD '12.
% J. Yang and J. Leskovec, Community-Affiliation Graph Model for 
% Overlapping Community Detection, ICDM '12.
% 
% The code works under Windows with Visual Studio or Cygwin with GCC,
% Mac OS X, Linux and other Unix variants with GCC. Make sure that a
% C++ compiler is installed on the system. Visual Studio project files
% and makefiles are provided. For makefiles, compile the code with
% "make all".
% 
% ////////////////////////////////////////////////////////////////////////
% Parameters:
%    -o: Output file name prefix (default:'')
%    -i: Input edgelist file name. DEMO: AGM with 2 communities
%    -l:  Input file name for node names (Node ID, Node label) 
%    -s: Random seed for AGM
%    -e: Edge probability between the nodes that do not share any 
%      community: set it to be 1 / N^2)
%    -c: Number of communities (0: determine it by AGM)
% 
% ////////////////////////////////////////////////////////////////////////
% Usage:
% 
% Detect 12 communities of universities (which correspond to NCAA 
% conferences) from the network of NCAA football teams:
% 
% agmfitmain -i:football.edgelist -l:football.labels -c:12 -e:0.1

ncommunities = 10 ;
% need to make labelfile
labnamefile = [OUTDIR_INTERM '/labnames_' PARC_STR '.txt' ] ;
writematrix([(1:NUM_NODES)' (1:NUM_NODES)'], ...
    labnamefile, ...
    'Delimiter','tab')

%% run svinet and recover communities

rng(123)
thr_vals = 0.05:0.01:0.15;

for thrIdx = 1:length(thr_vals) 

    dat_thr = thr_vals(thrIdx) ;

    % the file we want
    inputfile = [ OUTDIR_PROC '/groupavg_' PARC_STR ... 
            '_edgelist_thr' sprintf('%0.2f',dat_thr) '.txt' ] ;

    exe_call = strcat(binpath, ...
        ' -i:', inputfile, ...
        ' -l:', labnamefile, ...
        ' -c:', num2str(ncommunities) ) ;

    numRun = 100 ;
    agmfit_comms = cell(numRun,1) ;

    for idx = 1:numRun 
        disp(idx)

        tmpDir = strcat(OUTDIR_INTERM,'/tempDir/') ;
        mkdir(tmpDir) 
        cd(tmpDir)

        % run the command
        [~,~] = system([ exe_call ' -s:' num2str(randi(1000)) ' -o:' tmpDir ]) ;

        % communities data
        outCommPath = [ tmpDir '/cmtyvv.txt' ] ;
        outComm = dlmread(outCommPath) ;
        outComm = comms_mat_2_mat(outComm,NUM_NODES) ;
        agmfit_comms{idx} = outComm ; 

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
        write_comms(agmfit_comms{idx},[ tempdir '/comm' num2str(idx) '.txt'])
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

    agmfit_cent = agmfit_comms{centind} ;

%     %% look at it
% 
%     for idx = 1:numRun
%     imagesc(agmfit_comms{idx})
%     waitforbuttonpress
%     end

    %% make an agreement

    agree_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;
    agreew_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;

    for idx = 1:numRun
        disp(idx)

        agreew_dat(:,:,idx) = agmfit_comms{idx} * agmfit_comms{idx}' ;
        agree_dat(:,:,idx) = int32(agreew_dat(:,:,idx) > 0) ;

    end

    %% save results

    save([ OUTDIR_PROC '/agmfit_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ],...
        'agmfit_comms','agmfit_groups','agmfit_cent','centind') ;

end
