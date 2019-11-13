% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% load up some data

loadCons = load([ OUTDIR_PROC '/consensus_communities_all_' PARC_STR '.mat' ]) ;
loadNFC = load([ OUTDIR_PROC '/group+avg+fc+hcp-detrend.' PARC_STR '.rest.mat' ]) ;
loadEFC = load([ OUTDIR_PROC '/group_avg_hcp-detrend_rest_schaefer200-yeo17.mat' ]);
loadLabs = load([ OUTDIR_PROC '/hcp200.mat' ]) ;
loadEFCcomms = load([ OUTDIR_PROC '/consensus_communities_all_schaefer200-yeo17.mat' ]) ;

efcmat = loadEFC.rho + loadEFC.rho';
[~,sortLabs] = sort(loadLabs.lab) ;

nfcmat = loadNFC.rho + loadNFC.rho' ;

% get some indexing
n = NUM_NODES;
[ii,iii] = find(triu(ones(n),1));
einds = [ii,iii] ;

% colormap stuff
% cmapjet = jet(7);
% cmapjet(4,:) = 1;
% cmapjet = interp1(linspace(0,1,size(cmapjet,1)),cmapjet,linspace(0,1,256));

node_mask = logical(triu(ones(NUM_NODES),1)) ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% strengths
efc_str = strengths_und(efcmat)' ;
nefc_str = map_efcnodemeas_2_edges(efc_str,NUM_NODES) ;

% positive/negative strength
[efc_pos_str,efc_neg_str] = strengths_und_sign(efcmat) ;
nefc_pos_str = map_efcnodemeas_2_edges(efc_pos_str,NUM_NODES) ;
nefc_neg_str = map_efcnodemeas_2_edges(efc_neg_str,NUM_NODES) ;

% get the top/bottom 5% 
topefc = efc_pos_str>prctile(efc_pos_str, 95) ;
bttmefc = efc_neg_str>prctile(efc_neg_str, 95) ;
nefc_pos_top = map_efcnodemeas_2_edges(topefc,NUM_NODES) ;
nefc_neg_top = map_efcnodemeas_2_edges(bttmefc,NUM_NODES) ;

cmap = flipud(brewermap(256,'Spectral')) ;

% look at strength, look at pos/negative
figure
set(gcf, 'Units', 'Normalized', 'Position', [0, 0, 0.5, 0.5]);

ts = tight_subplot(2,3) ;
axes(ts(1))
imsc_grid_comm(nefc_str,loadLabs.lab,2,[ 0 0 0 ]) ; 
colormap(gca,cmap) 
colorbar; caxis([-250 250]), axis square ; set(gca,'YTickLabel',[]);
title('Overall strength')
axes(ts(2))
imsc_grid_comm(nefc_pos_str,loadLabs.lab,2,[ 1 1 1 ]) ;
colormap(gca,cmap(129:end,:)) 
colorbar; caxis([0 1500]); axis square ; set(gca,'YTickLabel',[]);
title('Positive strength')
axes(ts(3))
imsc_grid_comm(nefc_neg_str.*-1,loadLabs.lab,2,[ 1 1 1 ]) ; 
colormap(gca,cmap(1:128,:)) 
colorbar ; caxis([-1500 0]); axis square ; set(gca,'YTickLabel',[]);
title('Negative strength')

% relationship between node str and nefc str
axes(ts(4)) 
s = scatter(nfcmat(node_mask),nefc_str(node_mask),'filled') ;
s.MarkerFaceAlpha = 0.1 ; 
s.MarkerEdgeAlpha = 0 ;
xlabel('nFC edge weight')
ylabel('eFC overall strength')
axis square
cb = colorbar() ; cb.Visible = 'off' ;

cmap2 = [ 1 1 1 ; 0.2 0.2 0.2 ] ;

axes(ts(5))
imsc = imsc_grid_comm(nefc_pos_top+1,loadLabs.lab) ;
imsc.CDataMapping = 'direct' ;
colormap(gca,cmap2) 
caxis([0 1500]); axis square ; set(gca,'YTickLabel',[]);
cb = colorbar() ; cb.Visible = 'off' ;
title('Positive top %1')
axes(ts(6))
imsc = imsc_grid_comm(nefc_neg_top+1,loadLabs.lab) ;
imsc.CDataMapping = 'direct' ;
colormap(gca,cmap2)
caxis([0 1500]); axis square ; set(gca,'YTickLabel',[]);
title('Negative bottom %1')
cb = colorbar() ; cb.Visible = 'off' ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clustering coef
[efc_pos_cc,efc_neg_cc] = clustering_coef_wu_sign(efcmat) ;
nefc_pos_cc = map_efcnodemeas_2_edges(efc_pos_cc,NUM_NODES) ;
nefc_neg_cc = map_efcnodemeas_2_edges(efc_neg_cc,NUM_NODES) ;

% look at eeet
figure
ts = tight_subplot(1,2) ;
axes(ts(1))
imsc_grid_comm(nefc_pos_cc,loadLabs.lab) ; colorbar ; caxis([0 0.08]);
axis square ; set(gca,'YTickLabel',[]);
title('Positive clustering coef')
axes(ts(2))
imsc_grid_comm(nefc_neg_cc,loadLabs.lab) ; colorbar ; caxis([0 0.08]);
axis square ; set(gca,'YTickLabel',[]);
title('Negative clustering coef')

%% get system inds systems

% function [ comm_einds, within_inds, btwn_inds ] = efccomms_assignments(einds,ci)
[ yeo_einds , yeo_within , yeo_btwn ] = efccomms_assignments(einds,loadLabs.lab) ;

%% get some community metrics

% comm assignments for each edge, k=10, first layer is hcp
k10 = loadEFCcomms.Ci(:,10,1) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% participation coeff
[efc_pos_parti,efc_neg_parti] = participation_coef_sign(efcmat,k10) ;
nefc_pos_parti = map_efcnodemeas_2_edges(efc_pos_parti,NUM_NODES) ;
nefc_neg_parti = map_efcnodemeas_2_edges(efc_neg_parti,NUM_NODES) ;

% look at eeet
figure
ts = tight_subplot(1,2) ;
axes(ts(1))
imsc_grid_comm(nefc_pos_parti,loadLabs.lab) ; colorbar ; caxis([0.6 0.9]);
axis square ; set(gca,'YTickLabel',[]);
title('Positive participation coef')
axes(ts(2))
imsc_grid_comm(nefc_neg_parti,loadLabs.lab) ; colorbar ; caxis([0.6 0.9]);
axis square ; set(gca,'YTickLabel',[]);
title('Negative participation coef')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% module degree zscore
[efc_mdzs] = module_degree_zscore(efcmat,k10) ;
nefc_mdzs = map_efcnodemeas_2_edges(efc_mdzs,NUM_NODES) ;

% look at eeet
figure
ts = tight_subplot(1,2) ;
axes(ts(1))
imsc_grid_comm(nefc_mdzs,loadLabs.lab) ; colorbar ; caxis([-3 3]);
axis square ; set(gca,'YTickLabel',[]);
title('Module degree zscore')

% total participation
nefc_tot_parti = nefc_pos_parti + nefc_neg_parti ;

axes(ts(2))
s = scatter(nefc_tot_parti(node_mask),nefc_mdzs(node_mask),'filled') ;
s.MarkerFaceAlpha = 0.1 ; 
s.MarkerEdgeAlpha = 0 ;
axis square
xlabel('Overall participation coeff')
ylabel('Module degree zscore') 







