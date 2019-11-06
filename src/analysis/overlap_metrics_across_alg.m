% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% load up thresholded network 

dat_thr = 0.1 ;

loadDat = struct() ;

loadDat.nmf = load([ OUTDIR_PROC '/nmf_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
loadDat.svinet = load([ OUTDIR_PROC '/svinet_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
loadDat.agmfit = load([ OUTDIR_PROC '/agmfit_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
loadDat.thrlink = load([ OUTDIR_PROC '/thrlink_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;

methodstr = { 'nmf' 'svinet' 'agmfit' 'thrlink' } ;

% and load some data about the parcellation
lab_load = load([OUTDIR_PROC '/' PARC_STR_OTHER '.mat']) ;
labs_yeo = lab_load.lab ;

[~,sortyeo] = sort(labs_yeo) ;

%% load up previous entropy data

loadEnt = load([DROPBOX_DIR '/' DB_HCP_OUTPUT '/' PARC_STR '/' ENTROPY_RES ]) ;
loadCons = load([ OUTDIR_PROC '/consensus_communities_all_' PARC_STR '.mat' ]) ; 

%% gather dat

ovrRes = struct() ;
numRun = size(loadDat.agmfit.agmfit_comms,1) ;

for mmm = 1:length(methodstr)

    ovrRes.(methodstr{mmm}).agreew = zeros(NUM_NODES,NUM_NODES,numRun) ;
    ovrRes.(methodstr{mmm}).agree = zeros(NUM_NODES,NUM_NODES,numRun) ;
    ovrRes.(methodstr{mmm}).norment = zeros(NUM_NODES,numRun) ;
    
    % get comms variable name
    tmpStr = strcat(methodstr{mmm},'_comms') ;
    
    for idx = 1:numRun
        disp(idx)

        aa = loadDat.(methodstr{mmm}).(tmpStr){idx} * ...
            loadDat.(methodstr{mmm}).(tmpStr){idx}' ;
        ovrRes.(methodstr{mmm}).agreew(:,:,idx) = aa ;
        ovrRes.(methodstr{mmm}).agree(:,:,idx) = int32(aa > 0) ;

        % node entropies
        ovrRes.(methodstr{mmm}).norment(:,idx) = ...
            get_norm_node_ent(loadDat.(methodstr{mmm}).(tmpStr){idx}) ;
    end

    % make the mat for each
    ovrRes.(methodstr{mmm}).aw_mat = mean(ovrRes.(methodstr{mmm}).agreew,3) ;
    ovrRes.(methodstr{mmm}).a_mat = mean(ovrRes.(methodstr{mmm}).agree,3) ;
    
    % versatility from the agreement mat
    ovrRes.(methodstr{mmm}).vers = nodal_versatility(ovrRes.(methodstr{mmm}).a_mat) ;
    
    % average entropy
    ovrRes.(methodstr{mmm}).avg_ent = mean(ovrRes.(methodstr{mmm}).norment,2) ;
    
end

%% look at it

% organize by network
subplot(1,3,1)
imagesc(ovrRes.nmf.aw_mat(sortyeo,sortyeo))
subplot(1,3,2)
imagesc(ovrRes.agmfit.aw_mat(sortyeo,sortyeo))
subplot(1,3,3)
imagesc(ovrRes.svinet.aw_mat(sortyeo,sortyeo))


%%

efc_ent = mean(loadEnt.entnorm,2) ;
methodsCorr = corr([efc_ent ovrRes.thrlink.avg_ent ovrRes.svinet.avg_ent ovrRes.agmfit.avg_ent ovrRes.nmf.avg_ent]) ;







