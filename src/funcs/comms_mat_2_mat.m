function [ oMat ] = comms_mat_2_mat(iMat,nNodes) 

% read the number of communities
nComm = size(iMat,1) ;
% initialize the out
oMat = zeros(nNodes,nComm) ;

% loop over cell and put communities in columns
for idx = 1:nComm
    tmpInd = iMat(idx,:) ;
    oMat(tmpInd(tmpInd~=0),idx) = 1 ;
end


