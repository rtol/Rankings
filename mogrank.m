function [low mid high] = mogrank(cescore, bsscore, range, corr)
%[low mide high] = mogrank(cescore, bsscore)
%
%cescore is an n by m matrix of scores of n units for m bootstraps
%bsscore is the n by vector of actual scores
%range is the size of the confidence interval, in percentages
%
%mid returns the rank
%low and high return the bounds of the confidence interval specified by
%range
%
%the confidential interval of the rank is based on Mogstad, Romano,
%Shaikh and Wilhelm (2022 AEA P&P)
%
%if corr = false, the original procedure is followed
%if corr = true, the standard deviation of the difference is corrected for
%correlation
%
%Richard S.J. Tol, 20 May 2022


nobj = size(bsscore,1);
nbs = size(bsscore,2);
stderr = std(bsscore');
sd = stderr.^2 + stderr'.^2;
if corr
    covar = cov(bsscore');
    covar = covar - diag(diag(covar));
    sd = sd - 2*covar;
end
sd = sd.^0.5;
cediff = cescore' - cescore;

mcediff = zeros(nobj,nobj);
for i=1:nbs
    mcediff = mcediff + bsscore(:,i)'-bsscore(:,i);
end
mcediff = mcediff/nbs;

for i=1:nbs
    bsdiff = abs(cediff - (bsscore(:,i)'-bsscore(:,i)))./sd;
    cmat(:,i) = max(bsdiff)';
end

c = prctile(cmat',range)';

mid = rank(cescore);
low = 1 + sum(mcediff - sd.*repmat(c,[1 nobj]) > 0,2);
high = sum(mcediff + sd.*repmat(c,[1 nobj]) > 0,2);

end

