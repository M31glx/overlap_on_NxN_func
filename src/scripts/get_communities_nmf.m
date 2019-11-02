% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% load up thresholded network 

loadDat = load([ OUTDIR_PROC '/groupavg_' PARC_STR '_thr_networks.mat' ]) ;

thr_edgelists = loadDat.thr_edgelists ;
thr_networks = loadDat.thr_networks ;
thr_vals = loadDat.thr_vals ;

%% get comms

rng(123)

numRun = 100 ;
ncommunities = 10 ;

nmf_comms = cell(numRun,1) ;

inputMat = thr_networks{6} ;
dd = degrees_und(inputMat) ;
inputMat(~~eye(size(inputMat))) = dd ;

for idx = 1:numRun
    disp(idx)
    
    comm_bool = true ;
    while comm_bool
        % function [P,g,W,H] = commDetNMF(V,max_rank,W0,H0)
        [p,g] = commDetNMF(inputMat,ncommunities) ;
        
        %outComm = comms_cell_2_mat(g,NUM_NODES) ;
        outComm = p>0.25 ;

        if size(outComm,2) == ncommunities && (sum(sum(outComm) == 0) == 0)             
            comm_bool = false ;
        end
    end
 
    nmf_comms{idx} = outComm ;    
    
end

%% look at it

for idx = 1:numRun
imagesc(nmf_comms{idx})
waitforbuttonpress
end

 %% now get the centroid using 

% make a temp dir
tempdir = strcat(OUTDIR_INTERM,'/tempcomms/') ;
mkdir(tempdir)

% first write out all the files
for idx = 1:numRun
    disp(idx)
    write_comms(nmf_comms{idx},[ tempdir '/comm' num2str(idx) '.txt'])
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

nmf_cent = nmf_comms{centind} ;

%     %% look at it
% 
%     for idx = 1:numRun
%     imagesc(nmf_comms{idx})
%     waitforbuttonpress
%     end

%% make an agreement

agree_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;
agreew_dat = zeros(NUM_NODES,NUM_NODES,numRun) ;

for idx = 1:numRun
    disp(idx)

    agreew_dat(:,:,idx) = nmf_comms{idx} * nmf_comms{idx}' ;
    agree_dat(:,:,idx) = int32(agreew_dat(:,:,idx) > 0) ;

end

%% save results

save([ OUTDIR_PROC '/nmf_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ],...
    'nmf_comms','nmf_groups','nmf_cent','centind') ;


