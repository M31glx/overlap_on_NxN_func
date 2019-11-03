% clear some vars
clc ; clearvars ; close all
%% set the path

% run a config
config_file = 'config_schaefer200.m' ;
addpath(strcat(pwd,'/config'))
run(config_file) ;

%% read data

dat = load([ OUTDIR_PROC '/group+avg+fc+hcp-detrend.' PARC_STR '.rest.mat' ]) ;

% just ensure symmetric
rho = triu(dat.rho,1) + triu(dat.rho,1)' ;

%% get different thrs

thr_vals = 0.1:0.1:0.5;
thr_networks = cell(length(thr_vals),1) ;

% get min span tree
amst = graphminspantree(sparse(max(rho,[],'all') - rho),'method','kruskal') ;
amst = double(amst | amst');
amst = full(amst) ;

for idx = 1:length(thr_vals)

    % function c = fcn_mst_plus(a,amst,dens)
    thr_mask = fcn_mst_plus(rho,amst,thr_vals(idx)) ;
    thr_networks{idx} = rho .* thr_mask ;
end

%% convert to edge list
thr_edgelists = cell(length(thr_vals),1) ;

for idx = 1:length(thr_vals)
    
    net = thr_networks{idx} ; 
    net(net==0) = NaN ;
    el = Adj2Edg(net) ;
    thr_edgelists{idx} = el ;
end

%% write out the edge list

for idx = 1:length(thr_vals)
    disp(idx)
    filename = [ OUTDIR_PROC '/groupavg_' PARC_STR ... 
        '_edgelist_thr' sprintf('%0.2f',thr_vals(idx)) '.txt' ] ;
    writematrix(thr_edgelists{idx}(:,1:2),filename,'Delimiter','tab')
end

%% also save the thr networks to refer to in the future

save([ OUTDIR_PROC '/groupavg_' PARC_STR '_thr_networks.mat' ],...
    'thr_edgelists','thr_networks','thr_vals') 









