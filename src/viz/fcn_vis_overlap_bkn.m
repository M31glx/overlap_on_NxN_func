function ax = fcn_vis_overlap_bkn(theta,coor,sz,popout)

[n,k] = size(theta);
col = parula(k);
if nargin == 4
    ind = 1:k;
    ind(popout) = [];
    col(ind,:) = col(ind,:)*0.25;
end
thetanorm = bsxfun(@rdivide,theta,sum(theta,2));
thr = 0;
ax = axes;
hold(ax,'on');
for i = 1:n
    
    x = coor(i,1);
    y = coor(i,2);
    
    t = thetanorm(i,:);     % get the row
    t(t <= thr) = 0;        % floor any values below threshold
    t = t./sum(t);          % renormalize
    
    ind = find(t);          % get the indices of the nonzeros
    t = t(ind);             % and trim to just these inds
    m = length(t);          % get length of these inds
    
    if m == 1                                       % if only one community
        w = linspace(0,2*pi,100);
        a = cos(w)*sz(i) + x;
        b = sin(w)*sz(i) + y;
        fill(a,b,col(ind,:),'edgecolor','none');    % fill with color 
                                                    % according to comm
    else                                                
        tx = [0, cumsum(t)]*2*pi;
        for j = 1:m
            a = [0,cos(tx(j)),cos(linspace(tx(j),tx(j + 1),100)),cos(tx(j + 1)),0]*sz(i) + x;
            b = [0,sin(tx(j)),sin(linspace(tx(j),tx(j + 1),100)),sin(tx(j + 1)),0]*sz(i) + y;
            fill(a,b,col(ind(j),:),'edgecolor','none');
        end
    end
    
end
axis image off;