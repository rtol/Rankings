function [lowxie lower median upper upxie] = bsrankc(cescore, bsscore,range);

n = size(bsscore,1);
nbs = size(bsscore,2);

high = 100-(100-range)/2;
low = (100-range)/2;

stderr = std(bsscore');
sd = stderr.^2 + stderr'.^2;
covar = cov(bsscore');
sd = sd - 2*(covar - diag(diag(covar)));
stderr = sd.^0.5;

for i=1:n,
    rbs(i,:) = 0.5 + sum(normcdf(bsscore(:,:)-bsscore(i,:),0,repmat(stderr(i,:)',[1 nbs])),1);
end   

xie = normpdf(cescore-cescore',0,stderr);
xie = xie - diag(diag(xie));
xie(xie>1) = 1;
xie = sum(xie)';

median = prctile(rbs,50,2);
lower = prctile(rbs,low,2);
lowxie = floor(lower - 0.5*xie);
upper = prctile(rbs,high,2);
upxie = ceil(upper + 0.5*xie);

n = length(cescore);
lowxie(lowxie<1) = 1;
upxie(upxie>n) = n;

end