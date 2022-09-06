function ruy = crossCorrel(u,y,k)

n=length(u);

if size(u,1)>1 && size(u,2)==1
    u=u';
elseif size(u,1)>1
    error('u must be a vector.')
end

if size(y,1)>1 && size(y,2)==1
    y=y';
elseif size(u,1)>1
    error('u must be a vector.')
end

u_norm = u - mean(u);
y_norm = y - mean(y);

ruy = u_norm(1:n-k)*y_norm(1+k:n)'/(n-k);