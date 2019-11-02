function [] = write_comms(commMat,fname) 

nComm = size(commMat,2) ;

ff = fopen(fname,'w') ;

for idx = 1:nComm
    inds = find(commMat(:,idx))' ;
    fprintf(ff,'%i ',inds) ; fprintf(ff,'\n') ;
end

fclose(ff) ;
