function x = pop_decode(pop, pop_size, chromo_size, upper, lower)
x = zeros(pop_size, 1);
for i = 1:pop_size
    for j = 1:chromo_size
        if pop(i,j) == 1
            x(i) = x(i) + 2^(j-1);
        end
    end
    x(i) = (x(i) / (2^chromo_size-1)) * (upper - lower) + lower;
end