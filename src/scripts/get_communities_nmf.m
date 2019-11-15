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

ncommunities = 10 ;
nmf_comms = cell(NUM_RUN,1) ;

rng(123)
%thr_vals = 0.05:0.01:0.15;
%thr_vals = 0.1 ;

for thrIdx = 1:length(thr_vals) 

    dat_thr = thr_vals(thrIdx) ;
    
    inputMat = thr_networks{thrIdx} ;
    dd = degrees_und(inputMat) ;
    inputMat(~~eye(size(inputMat))) = dd ;

    for idx = 1:NUM_RUN
        disp(idx)

        comm_bool = true ;
        while comm_bool
            % function [P,g,W,H] = commDetNMF(V,max_rank,W0,H0)
            [p,g] = commDetNMF(inputMat,ncommunities) ;

            %outComm = comms_cell_2_mat(g,NUM_NODES) ;
            outComm = p>COMM_PROB_THR ;

            if size(outComm,2) == ncommunities && (sum(sum(outComm) == 0) == 0)             
                comm_bool = false ;
            end
        end

        nmf_comms{idx} = outComm ;    

    end

    %% look at it

    % for idx = 1:NUM_RUN
    % imagesc(nmf_comms{idx})
    % waitforbuttonpress
    % end

     %% now get the centroid using 

    % now compare all communities
    ovrmutinfomat = zeros(NUM_RUN) ;
    for idx = 1:NUM_RUN
       for jdx = 1:NUM_RUN
           if idx >= jdx
               continue
           else
               disp([num2str(idx) ' ' num2str(jdx)])
           end
           
           ovrmutinfomat(idx,jdx) = gnmi(nmf_comms{idx},nmf_comms{jdx},NUM_NODES) ;
           
        end
    end

    %  find the centroid
    nmf_simmat = ovrmutinfomat + ovrmutinfomat' ;
    [~,centind] = max(sum(nmf_simmat)) ;

    nmf_cent = nmf_comms{centind} ;

    %     %% look at it
    % 
    %     for idx = 1:NUM_RUN
    %     imagesc(nmf_comms{idx})
    %     waitforbuttonpress
    %     end

    %% make an agreement

    agree_dat = zeros(NUM_NODES,NUM_NODES,NUM_RUN) ;
    agreew_dat = zeros(NUM_NODES,NUM_NODES,NUM_RUN) ;

    for idx = 1:NUM_RUN
        disp(idx)

        agreew_dat(:,:,idx) = nmf_comms{idx} * nmf_comms{idx}' ;
        agree_dat(:,:,idx) = int32(agreew_dat(:,:,idx) > 0) ;

    end

    %% save results

    save([ OUTDIR_PROC '/nmf_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ],...
        'nmf_*','centind') ;

end


