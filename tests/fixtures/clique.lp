obj: x[1] + x[2] + x[3] + x[4]
subject to
 clq1: x[1] + x[2] + x[3] <= 1
 clq2: x[3] + x[4] <= 1
bounds
 0 <= x[1] <= 1
 0 <= x[2] <= 1
 0 <= x[3] <= 1
 0 <= x[4] <= 1
end
