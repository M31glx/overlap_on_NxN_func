% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% load up thresholded network 

thresholds = [ 0.1 0.2 0.3 0.4 ] ;
loadDat = struct() ;

% methodstr = { 'nmf' 'svinet' 'agmfit' 'thrlinkbin' 'thrlink' } ;
% 'thrlinkbin' is effectively the same as 'thrlink' method
methodstr = { 'nmf' 'svinet' 'agmfit' 'thrlink' } ;

for idx = 1:length(thresholds) 

    disp(idx)
    dat_thr = thresholds(idx) ;

    for jdx = 1:length(methodstr)
        disp(jdx)

        loadDat(idx).(methodstr{jdx}) = load([ OUTDIR_PROC '/' ...
            methodstr{jdx} '_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ]) ;
    end   
end

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
        ovrRes(tdx).(methodstr{mmm}).numcomms = zeros(NUM_NODES,numRun) ;
        
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
            
            % number comms
            ovrRes(tdx).(methodstr{mmm}).numcomms(:,idx) = ...
                sum(loadDat(tdx).(methodstr{mmm}).(tmpStr){idx},2) ;
            
        end

        % make the mat for each
        ovrRes(tdx).(methodstr{mmm}).aw_mat = mean(ovrRes(tdx).(methodstr{mmm}).agreew,3) ;
        ovrRes(tdx).(methodstr{mmm}).a_mat = mean(ovrRes(tdx).(methodstr{mmm}).agree,3) ;

        % versatility from the agreement mat
        ovrRes(tdx).(methodstr{mmm}).vers = nodal_versatility(ovrRes(tdx).(methodstr{mmm}).a_mat) ;

        % average entropy
        ovrRes(tdx).(methodstr{mmm}).avg_ent = mean(ovrRes(tdx).(methodstr{mmm}).norment,2) ;

        % number of communities per method
        ovrRes(tdx).(methodstr{mmm}).avg_numc = median(ovrRes(tdx).(methodstr{mmm}).numcomms,2) ;
            
    end

end

efc_ent = mean(loadEnt.entnorm,2) ;
efc_numc = median(loadEnt.tot,2) ;

%%

plotit = true ;

capitl = @(x) [upper(x(1)) x(2:end)] ;

% function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
if ~plotit
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.4 0.6]);      
ts = tight_subplot(1,length(thresholds)) ;
end

for idx = 1:length(thresholds)

    if ~plotit
    axes(ts(idx))
    end
    
    methodsCorr = corr([efc_ent ...
        ovrRes(idx).thrlink.avg_ent ...
        ovrRes(idx).agmfit.avg_ent ...
        ovrRes(idx).svinet.avg_ent ...
        ovrRes(idx).nmf.avg_ent]) ;
    imsc = imagesc(methodsCorr) ;

    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel', [ 'eFC' fliplr(cellfun(capitl,methodstr,'UniformOutput',false))])
    xtickangle(45)
    colorbar()
    caxis([0 1])
    title(['Threshold ' num2str(thresholds(idx))])
    
    set(gca,'ytick',[])
    axis square

    if plotit
        print(gcf,...
        '-dpdf',[ PDF_FIG_DIR '/thr_' num2str(thresholds(idx)) '_similaritymat.pdf' ]);
        close(gcf);
    end
    
end

%%

plotit = true ;

capitl = @(x) [upper(x(1)) x(2:end)] ;

rng(123)

if ~plotit
% function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5]);      
ts = tight_subplot(1,length(thresholds)) ;
end


for kdx = 1:length(thresholds)
for idx = fliplr(1:4)

    thridx = kdx ;
    
    if ~plotit
    axes(ts(5-idx))
    end
    
    scatcolros = loadCmap.cmap(loadAff.lab) ;
    colormap(loadCmap.cmap)

    ss = scatter(ovrRes(thridx).(methodstr{idx}).avg_ent,...
                 efc_ent,...
                 [],scatcolros,...
                 'filled','MarkerEdgeColor',[.1 .1 .1]) ;
    ss.MarkerFaceAlpha = 0.8 ;
    ss.MarkerEdgeAlpha = 0.1 ;
    ylim([-0.05 1.05])
    xlim([-0.05 1.05])
    axis square

    if ~plotit
    if idx == 4
    ylabel('eFC entropy')
    else
       set(gca,'YTickLabel',[]);
    end
    else
       ylabel('eFC entropy') 
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
    text(xr(2)*0.7,.1,annotText,'fontsize',14,'Interpreter','latex')

    if plotit
        print(gcf,...
        '-dpdf',[ PDF_FIG_DIR '/thr_' num2str(thresholds(kdx)) '_scatter_efc_vs_' methodstr{idx} '.pdf' ]);
        close(gcf);
    end

end
end

%% histogram of median number of comms per node
% 
% 
%
% figure 
% edges = (1:15)-0.5 ;
% 
% %ts = tight_subplot(1,length(thresholds)) ;
% subplot(1,length(thresholds),1) ;
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5 0.5]);      
% 
% for idx = 1:length(thresholds)
% 
%     subplot(1,length(thresholds),idx) 
% 
%     hold on
%     
%     for jdx = 1:length(methodstr)
% 
%         histogram(ovrRes(idx).(methodstr{jdx}).avg_numc,edges)
% 
%     end
% 
%     histogram(efc_numc,edges)
% %     
%     hold off
% 
% end

%% package data

overlapMetricOthrAlgs = ovrRes ;
overlapComms = loadDat ;
efcEnt = loadEnt ;

save([OUTDIR_PROC '/overlap_metrics_otherAlgs_' PARC_STR '.mat'],...
    'overlapMetricOthrAlgs','overlapComms',...
    'thresholds','efcEnt','-v7.3')

%% super simiple package data

overlap_norment = struct() ;

for idx = 1:length(thresholds)
    
    overlap_norment(idx).eFC = efc_ent ;
    
    for jdx = 1:length(methodstr)
        disp(jdx)
        overlap_norment(idx).(methodstr{jdx}) = ovrRes(idx).(methodstr{jdx}).avg_ent ;
    end

end

save([OUTDIR_PROC '/overlap_norment' PARC_STR '.mat'],...
    'overlap_norment', 'thresholds','-v7.3')


