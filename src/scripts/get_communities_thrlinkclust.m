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

%% inline func

% https://www.mathworks.com/help/matlab/math/floating-point-numbers-within-specific-range.html
randsamprange = @(minn,maxx,numm)(maxx-minn).*rand(numm,1) + minn;

%% get comms

rng(123)

ncommunities = 10 ;
thrlink_comms = cell(NUM_RUN,1) ;
thrlink_wcomms = cell(NUM_RUN,1) ;

%thr_vals = 0.05:0.01:0.15;
%thr_vals = 0.1 ;

for thrIdx = 2:length(thr_vals) 
    
    % make input mat into a line graph, and then make into mod mat
    [lg,~,u,v] = fcn_line_graph(thr_networks{thrIdx}); % calculate line graph
    m = size(lg,1) ;

    % first need to do gamma sweep
    initgamas = 0.1:0.01:6 ;
    sweepComs = zeros(m,length(initgamas)) ; 
    for idx = 1:length(initgamas)
         % newman girvan mod
        b = modularity(lg,initgamas(idx)) ;       
        sweepComs(:,idx) = genlouvain(b,[],false);
        % early stopping criteria
        tmpnumcoms = max(sweepComs(:,idx)) ;
        disp(tmpnumcoms)
        if tmpnumcoms >= ncommunities + 5
            break
        end
    end

    % get num coms at each sweep iteration 
    numCrossGam = max(sweepComs)' ;
    kInd = find(numCrossGam == ncommunities) ;

    if isempty(kInd) || (kInd(1) == kInd(end))
       error('gamma range not big enough') 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % now run it 
    
    % new gamma range
    startG = initgamas(kInd(1)) ; endG = initgamas(kInd(end)) ;
    gammas = randsamprange(startG,endG,NUM_RUN) ;
    
    % preallocations
    ci = zeros(m,NUM_RUN);            % preallocate for communities
    q = zeros(NUM_RUN,1);             % preallocate for modularity, q
    ent = zeros(NUM_NODES,NUM_RUN); % preallocate for entropy
    entnorm = ent;                  % normalized entropy
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    idx=1 ;
    lattempt = 1 ;
    maxlattempt = 5 ; 
    while idx <= NUM_RUN
        disp(idx)
               
        if lattempt > maxlattempt
            % get a new gamma here
            gammas(idx) = randsamprange(startG,endG,1) ;
            lattempt = 1 ;
        end
        
        b = modularity(lg,gammas(idx)) ;       
        [ci(:,idx),q(idx)] = genlouvain(b,[],false);
        
        if ci(:,idx) ~= ncommunities
           lattempt = lattempt + 1 ;
           continue
        else % sucessfull iter
           disp([ 'generated com:' num2str(idx) ]) 
        end
        
        [ent(:,idx),entnorm(:,idx),h] = ...   % calculate entropies and comms
            fcn_node_entropy(ci(:,idx),u,v,NUM_NODES);
                
        thrlink_wcomms{idx} = h ;
        thrlink_comms{idx} = h > 0 ;    

        % increment
        idx = idx + 1 ;
    end

    %% look at it

    % for idx = 1:NUM_RUN
    % imagesc(thrlink_comms{idx})
    % waitforbuttonpress
    % end

     %% now get the centroid using 

    % make a temp dir
    tempdir = strcat(OUTDIR_INTERM,'/tempcomms/') ;
    mkdir(tempdir)

    % first write out all the files
    for idx = 1:NUM_RUN
        disp(idx)
        write_comms(thrlink_comms{idx},[ tempdir '/comm' num2str(idx) '.txt'])
    end

    % now compare all communities
    ovrmutinfomat = zeros(NUM_RUN) ;
    mutexe = [ PROJECT_DIR '/src/external/mutual3/mutual' ] ;
    for idx = 1:NUM_RUN
       for jdx = 1:NUM_RUN
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
    thrlink_simmat = ovrmutinfomat + ovrmutinfomat' ;
    [~,centind] = max(sum(thrlink_simmat)) ;

    thrlink_cent = thrlink_comms{centind} ;

    %     %% look at it
    % 
    %     for idx = 1:NUM_RUN
    %     imagesc(thrlink_comms{idx})
    %     waitforbuttonpress
    %     end

    %% make an agreement

    agree_dat = zeros(NUM_NODES,NUM_NODES,NUM_RUN) ;
    agreew_dat = zeros(NUM_NODES,NUM_NODES,NUM_RUN) ;

    for idx = 1:NUM_RUN
        disp(idx)

        agreew_dat(:,:,idx) = thrlink_comms{idx} * thrlink_comms{idx}' ;
        agree_dat(:,:,idx) = int32(agreew_dat(:,:,idx) > 0) ;

    end

    %% save results

    dat_thr = thr_vals(thrIdx) ;
    
    save([ OUTDIR_PROC '/thrlink_overlapcomms_' PARC_STR '_thr' sprintf('%0.2f',dat_thr) '_networks.mat' ],...
        'thrlink_*','centind') ;

end

