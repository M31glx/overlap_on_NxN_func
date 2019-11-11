function [ comm_einds, within_inds, btwn_inds ] = efccomms_assignments(einds,ci)

% community einds
comm_einds = [ci(einds(:,1)) ci(einds(:,2)) ] ;

% when edge compared to edge in same community 
within_inds = zeros(length(einds),max(ci)) ;
% given a community, edges connecting to community (not within)
btwn_inds = zeros(length(einds),max(ci)) ;

% loop through communities
for idx = 1:max(ci)
    within_inds(:,idx) = comm_einds(:,1)==idx & comm_einds(:,2) == idx ;
    btwn_inds(:,idx) = (comm_einds(:,1)==idx & comm_einds(:,2) ~= idx) | ...
        (comm_einds(:,2)==idx & comm_einds(:,1) ~= idx) ;    
end
