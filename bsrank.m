function [lowxie low mid high highxie] = bsrank(cescore, bsscore,range);
%[lowxie lower rank upper upxie] = bsrank(score,range,bw,f);
%
%cescore is an n by 1 vector of scores
%bsscore is an n by m matrix of scores of n units for m bootstraps
%range is the size of the confidence interval, in percentages
%
%mid returns the rank of score
%
%low and high return the bootstrap bounds of the confidence interval,
%following Goldstein and Spiegelhalter (1996 JRSS-A)
%https://rss.onlinelibrary.wiley.com/doi/10.2307/2983325
%
%lowxie and highxie return the bootstrap bounds corrected for near-ties,
%following Xie, Singh and Zhang (2012 JASA)
%https://www.tandfonline.com/doi/abs/10.1198/jasa.2009.0142
%
%Richard S.J. Tol, 20 May 2022

if range > 100 | range < 0
    disp('confidence interval set to 95%')
    range = 95;
end

up = 100-(100-range)/2;
down = (100-range)/2;

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
    
mid = rank(cescore);
low = prctile(rbs,down,2);
lowsm = prctile(rbsx,down,2);
lowxie = floor(lowsm - 0.5*xie);
high = prctile(rbs,up,2);
highsm = prctile(rbsx,up,2);
highxie = ceil(highsm + 0.5*xie);

n = length(cescore);
lowxie(lowxie<1) = 1;
highxie(highxie>n) = n;

end