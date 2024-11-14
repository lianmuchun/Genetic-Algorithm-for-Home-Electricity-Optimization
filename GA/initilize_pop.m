function pop = initilize_pop(pop_size, chromo_size)
pop = zeros(pop_size, chromo_size);
for i = 1:pop_size
    for j = 1:chromo_size
        pop(i,j) = round(rand);
    end
end