\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[7] + x[3] + x[8] + x[4] + x[5] + x[6]
Subject To
 clq1: x[1] + x[2] + x[3] <= 1
 clq2: x[7] + x[5] + x[6] <= 1
 clq3: x[1] + x[7] + x[8] <= 1
 clq4: x[2] + x[4] <= 1
 clq5: x[3] + x[4] + x[5] <= 1
 clq6: x[8] + x[6] <= 1
 clq7: x[2] + x[8] <= 1
 clq8: x[4] + x[6] <= 1
Bounds
Binaries
 x[1] x[2] x[7] x[3] x[8] x[4] x[5] x[6]
End
