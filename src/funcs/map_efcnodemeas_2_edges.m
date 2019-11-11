function [ outMat ] = map_efcnodemeas_2_edges(efcnodemeas,numnodes)

% get inds
[ii,iii] = find(triu(ones(numnodes),1));

% initialize out mat
outMat = zeros(numnodes,numnodes);
% make self loops NaN
outMat(~~eye(numnodes)) = NaN ;

% for each node
for idx = 1:numnodes
    % get the N-1 edges
    ind = ii == idx | iii == idx ; 
    % write out the row
    outind = ~isnan(outMat(idx,:)) ;
    % the N-1 into the outMat
    outMat(idx,outind) = efcnodemeas(ind) ;
end
