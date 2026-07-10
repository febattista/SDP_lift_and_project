obj: x[1] + x[2] + x[3]
subject to
 nod1: x[2] + x[3] + 2 x[1] <= 2
 edge1: x[2] + x[3] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
 0 <= x[3] <= 1
end
