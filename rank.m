function rbs = rank(score)

for i=1:319,
    rbs(i) = 0.5 + sum(score(i)<score) + 0.5*sum(score(i)==score);
end

end

