function [c median lower upper] = mogrank(cescore, bsscore, range, corr)
%[c rank lower upper] = mogrank(cescore, bsscore)
%
%cescore is an nxm matrix of scores of n units for m bootstraps
%bsscore is the nx1 vector of actual scores
%range is the size of the confidence interval, in percentages
%
%median returns the median rank
%lower and upper return the bounds of the confidence interval specified by
%range
%
%the confidential interval of the rank is based on Mogstad, Romano, Shaikh and Wilhelm (2022
%AEA P&P)
%
%if corr = true, the original procedure is followed
%if corr = false, the standard deviation of the difference is corrected for
%correlation

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

median = 1 + sum(mcediff > 0,2);
lower = 1 + sum(mcediff - sd.*repmat(c,[1 nobj]) > 0,2);
upper = sum(mcediff + sd.*repmat(c,[1 nobj]) > 0,2);

end

