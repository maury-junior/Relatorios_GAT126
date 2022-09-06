function ruu = autoCorrel(u,k)

n=length(u);

if size(u,1)>1 && size(u,2)==1
    u=u';
elseif size(u,1)>1
    error('u must be a vector.')
end

u_norm = u - mean(u);

num = u_norm(1:n-k)*u_norm(1+k:n)';
den = u_norm*u_norm';

ruu = num/den;