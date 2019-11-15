% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% load up thresholded network 

thresholds = [ 0.1 0.2 0.3 0.4 0.5 ] ;
loadDat = struct() ;
    
for idx = 1:length(thresholds) 

    disp(idx)
    dat_thr = thresholds(idx) ;

    loadDat(idx).nmf = load([ OUTDIR_PROC '/nmf_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
    loadDat(idx).svinet = load([ OUTDIR_PROC '/svinet_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
    loadDat(idx).agmfit = load([ OUTDIR_PROC '/agmfit_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
    loadDat(idx).thrlink = load([ OUTDIR_PROC '/thrlink_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;

end

methodstr = { 'nmf' 'svinet' 'agmfit' 'thrlink' } ;

% and load some data about the parcellation
lab_load = load([OUTDIR_PROC '/' PARC_STR_OTHER '.mat']) ;
labs_yeo = lab_load.lab ;

[~,sortyeo] = sort(labs_yeo) ;

%% load up previous entropy data

loadEnt = load([DROPBOX_DIR '/' DB_HCP_OUTPUT '/' PARC_STR '/' ENTROPY_RES ]) ;
loadCons = load([ OUTDIR_PROC '/consensus_communities_all_' PARC_STR '.mat' ]) ; 
loadCmap = load([ PROJECT_DIR '/data/processed/hcp200cmapOther.mat' ]) ;
loadAff = load([ PROJECT_DIR '/data/processed/hcp200.mat' ]) ;

%% gather dat

ovrRes = struct() ;
numRun = size(loadDat(1).agmfit.agmfit_comms,1) ;

for tdx = 1:length(thresholds)

    for mmm = 1:length(methodstr)

        ovrRes(tdx).(methodstr{mmm}).agreew = zeros(NUM_NODES,NUM_NODES,numRun) ;
        ovrRes(tdx).(methodstr{mmm}).agree = zeros(NUM_NODES,NUM_NODES,numRun) ;
        ovrRes(tdx).(methodstr{mmm}).norment = zeros(NUM_NODES,numRun) ;

        % get comms variable name
        tmpStr = strcat(methodstr{mmm},'_comms') ;

        for idx = 1:numRun
            disp(idx)

            aa = loadDat(tdx).(methodstr{mmm}).(tmpStr){idx} * ...
                loadDat(tdx).(methodstr{mmm}).(tmpStr){idx}' ;
            ovrRes(tdx).(methodstr{mmm}).agreew(:,:,idx) = aa ;
            ovrRes(tdx).(methodstr{mmm}).agree(:,:,idx) = int32(aa > 0) ;

            % node entropies
            ovrRes(tdx).(methodstr{mmm}).norment(:,idx) = ...
                get_norm_node_ent(loadDat(tdx).(methodstr{mmm}).(tmpStr){idx}) ;
        end

        % make the mat for each
        ovrRes(tdx).(methodstr{mmm}).aw_mat = mean(ovrRes(tdx).(methodstr{mmm}).agreew,3) ;
        ovrRes(tdx).(methodstr{mmm}).a_mat = mean(ovrRes(tdx).(methodstr{mmm}).agree,3) ;

        % versatility from the agreement mat
        ovrRes(tdx).(methodstr{mmm}).vers = nodal_versatility(ovrRes(tdx).(methodstr{mmm}).a_mat) ;

        % average entropy
        ovrRes(tdx).(methodstr{mmm}).avg_ent = mean(ovrRes(tdx).(methodstr{mmm}).norment,2) ;

    end

end

efc_ent = mean(loadEnt.entnorm,2) ;

%% look at it

% % organize by network
% subplot(1,4,1)
% imagesc(ovrRes.nmf.aw_mat(sortyeo,sortyeo))
% subplot(1,4,2)
% imagesc(ovrRes.agmfit.aw_mat(sortyeo,sortyeo))
% subplot(1,4,3)
% imagesc(ovrRes.svinet.aw_mat(sortyeo,sortyeo))
% subplot(1,4,4)
% imagesc(ovrRes.thrlink.aw_mat(sortyeo,sortyeo))
% 
% subplot(1,4,1)
% imagesc(ovrRes.nmf.a_mat(sortyeo,sortyeo))
% subplot(1,4,2)
% imagesc(ovrRes.agmfit.a_mat(sortyeo,sortyeo))
% subplot(1,4,3)
% imagesc(ovrRes.svinet.a_mat(sortyeo,sortyeo))
% subplot(1,4,4)
% imagesc(ovrRes.thrlink.a_mat(sortyeo,sortyeo))

%%

capitl = @(x) [upper(x(1)) x(2:end)] ;

% function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.6]);      
ts = tight_subplot(1,length(thresholds)) ;

for idx = 1:length(thresholds)

    axes(ts(idx))
    
    methodsCorr = corr([efc_ent ...
        ovrRes(idx).thrlink.avg_ent ovrRes(idx).agmfit.avg_ent ...
        ovrRes(idx).svinet.avg_ent ovrRes(idx).nmf.avg_ent]) ;
    imsc = imagesc(methodsCorr) ;

    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel', [ 'eFC' fliplr(cellfun(capitl,methodstr,'UniformOutput',false))])
    xtickangle(45)
    colorbar()
    caxis([-0.5 1])
    title(['Threshold ' num2str(thresholds(idx))])

    axis square

end

%%

capitl = @(x) [upper(x(1)) x(2:end)] ;

rng(123)

% function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5]);      
ts = tight_subplot(1,4) ;

thridx = 1 ;

for idx = fliplr(1:4)
    
    axes(ts(5-idx))

    scatcolros = loadCmap.cmap(loadAff.lab) ;
    colormap(loadCmap.cmap)

    ss = scatter(ovrRes(thridx).(methodstr{idx}).avg_ent,...
                 efc_ent,...
                 [],scatcolros,...
                 'filled','MarkerEdgeColor',[.1 .1 .1]) ;
    ss.MarkerFaceAlpha = 0.8 ;
    ss.MarkerEdgeAlpha = 0.1 ;
    ylim([0 1])
    xlim([0 1])
    axis square

    if idx == 4
    ylabel('eFC entropy')
    else
       set(gca,'YTickLabel',[]);
    end
    xlabel([ capitl(methodstr{idx}) ' entropy'])
    
    cc= corr(ovrRes(thridx).(methodstr{idx}).avg_ent,efc_ent,'type','spearman') ;

%     % quick bootstrap?
%     nperm = 500 ; 
%     randinds = randi(NUM_NODES,[ NUM_NODES nperm]) ;
%     bootcc = zeros(nperm,1) ;
%     for bdx = 1:nperm
%         bootcc(bdx) = corr(ovrRes.(methodstr{idx}).avg_ent(randinds(:,bdx)),efc_ent(randinds(:,bdx)),'type','spearman') ;
%     end
%     ci = prctile(bootcc,[2.5 97.5]) ;
    
    annotText = [ '$\rho = ' num2str(round(cc,2)) '$' ] ;
    % '[' num2str(round(ci(1),2)) ',' num2str(round(ci(2),2)) ']'
    xr = xlim ; yr = ylim ;
    text(0.6,.1,annotText,'fontsize',14,'Interpreter','latex')

    
end

%% package data

overlapMetricOthrAlgs = ovrRes ;
overlapComms = loadDat ;
efcEnt = loadEnt ;

save([OUTDIR_PROC '/overlap_metrics_otherAlgs_' PARC_STR '.mat'],...
    'overlapMetricOthrAlgs','overlapComms',...
    'thresholds','efcEnt','-v7.3')

