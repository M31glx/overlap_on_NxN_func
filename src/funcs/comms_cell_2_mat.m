function [ oMat ] = comms_cell_2_mat(iCell,nNodes) 

% read the number of communities
nComm = length(iCell) ;
% initialize the out
oMat = zeros(nNodes,nComm) ;

% loop over cell and put communities in columns
for idx = 1:nComm
    tmpInd = iCell{idx} ;
    oMat(tmpInd,idx) = 1 ;
end


