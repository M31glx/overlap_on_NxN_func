function [ gammarange ] = sweep_gamma(W,k,initgamas,modfunc,kbuff)

if nargin < 2
   error('need at least 2 gammas') 
end

if ~exist('initgamas','var') || isempty(initgamas)
    initgamas = 0.1:0.01:10 ;
end

if ~exist('kbuff','var') || isempty(kbuff)
    kbuff = 4 ;
end

if ~exist('modfunc','var') || ~isa(modfunc,'function_handle') || isempty(modfunc)
    disp('did not recognize a modularity func, using newman-girvan')
    modfunc = @modularity_ng ;
else
    disp('using provided modularity func')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first need to do gamma sweep
sweepComs = zeros(size(W,1),length(initgamas)) ; 

for idx = 1:length(initgamas)
     % newman girvan mod
    b = modfunc(W,initgamas(idx)) ;       
    sweepComs(:,idx) = genlouvain(b,[],false);
    % early stopping criteria
    tmpnumcoms = max(sweepComs(:,idx)) ;
    disp(tmpnumcoms)
    if tmpnumcoms >= k + kbuff
        break
    end
end

% get num coms at each sweep iteration 
numCrossGam = max(sweepComs)' ;
kInd = find(numCrossGam == k) ;

if isempty(kInd) || (kInd(1) == kInd(end))
   error('gamma range not big enough') 
end

gammarange = [ initgamas(kInd(1)) initgamas(kInd(end)) ] ;

end

function [B] = modularity_ng(A,gamma)
%MODULARITY returns monolayer Newman-Girvan modularity matrix for network given by adjacency matrix A, matrix version
%
% Version: 2.2.0
% Date: Thu 11 Jul 2019 12:25:42 CEST
%
if nargin<2||isempty(gamma)
	gamma=1;
end
k=sum(A,2);
d=sum(A,1);
twom=sum(k);
B=full((A+A')/2-gamma/2*(k*d+d'*k')/twom);
end

