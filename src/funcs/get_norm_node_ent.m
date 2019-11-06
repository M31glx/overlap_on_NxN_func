function [ normNodeEnt ] = get_norm_node_ent(iComs) 

probs = bsxfun(@rdivide,iComs,sum(iComs,2));
nodeEnt = -nansum(probs.*log2(probs),2);
normNodeEnt = nodeEnt/log2(size(iComs,2)) ;
