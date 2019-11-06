function [ nodal_vers ] = nodal_versatility(iAgreeMat)

% check to make sure input mat is already scaled
if min(iAgreeMat) < 0 | max(iAgreeMat) > 1
   error('input mat needs to be scaled between 0 and 1') 
end

sineAgreeMat = sin(pi*iAgreeMat);
nodal_vers = sum(sineAgreeMat, 1)./sum(iAgreeMat, 1);
nodal_vers(nodal_vers<1e-10) = 0;
nodal_vers = nodal_vers(:) ;