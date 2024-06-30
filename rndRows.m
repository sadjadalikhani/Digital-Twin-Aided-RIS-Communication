function [r] = rndRows(grid,grid_num)
    rng("shuffle")
    r = floor(grid(1) + (grid(2)-grid(1)) .* rand(1,grid_num));
end