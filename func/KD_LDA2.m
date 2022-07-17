function v=KD_LDA2(rr,ll,rs)

 
if nargin==3;
    x = [ll;rr;rs];
    groupVar = [ones(size(ll,1),1); ones(size(rr,1),1)*2; ones(size(rs,1),1)*3];
else
    x = [ll;rr];
    groupVar = [ones(size(ll,1),1); ones(size(rr,1),1)*2];
end

xm = mean(x);
n=size(x,1);
x = x - xm(ones(n,1),:);       % done with the original x
T = x'*x;

% Now compute the Within sum of squares matrix
W = zeros(size(T));
for j=1:2
    r = find(groupVar == j);
    nr = length(r);
    if (nr > 1)
        z = x(r,:);
        xm = mean(z);
        z = z - xm(ones(nr,1),:);
        W = W + z'*z;
    end
end

B=T-W;

 
[u,s,v]=svd(pinv(W)*B);

cr=rr*v;
cl=ll*v;
sgn=sign(mean(cr(:,1))-mean(cl(:,1)));
v(:,1)=v(:,1)*sgn;
