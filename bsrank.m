function [lowxie lower median upper upxie] = bsrank(cescore, bsscore,range);
%[lower median upper] = bsrank(score,range,bw,f);
%
%score is an nxm matrix of scores of n units for m bootstraps
%range is the size of the confidence interval, in percentages
%bw is the bandwidth
%corr is an optional parameter
%
%median returns the median rank
%lower and upper return the bounds of the confidence interval specified by
%range
%
%if bw = 0, the confidential interval of the rank is based on Goldstein and Spiegelhalter (1996
%JRSS-A) https://rss.onlinelibrary.wiley.com/doi/10.2307/2983325
%
%if bw > 0, the confidential interval is corrected for ties following Xie,
%Singh and Zhang (2012 JASA)
%https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.0142 where the
%bandwidth equals bw times the bootstrap standard deviation of the ranked
%object

if range > 100 | range < 0
    disp('confidence interval set to 95%')
    range = 95;
end

high = 100-(100-range)/2;
low = (100-range)/2;

iqr = prctile(cescore,75)-prctile(cescore,25);

n = size(bsscore,1);
stderr = std(bsscore');

for i=1:n,
    rbs(i,:) = 0.5 + sum(bsscore(:,:)>bsscore(i,:),1) + 0.5*sum(bsscore(:,:)==bsscore(i,:),1);
    rbsx(i,:) = 0.5 + sum(normcdf(bsscore(:,:)-bsscore(i,:),0,iqr*sqrt(stderr(i))),1);
end

xie = normpdf(cescore-cescore',0,iqr*sqrt(stderr)');
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