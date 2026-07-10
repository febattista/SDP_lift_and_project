obj: x[1] + x[2] + x[3]
subject to
 edge1: x[1] + x[2] <= 1
 edge2: x[1] + x[3] <= 1
 edge3: x[2] + x[3] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
 0 <= x[3] <= 1
end
