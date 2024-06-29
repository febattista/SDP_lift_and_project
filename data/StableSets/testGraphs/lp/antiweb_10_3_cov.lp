\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[9] + x[3] + x[10] + x[4] + x[5] + x[6] + x[7] + x[8]
Subject To
 clq1: x[1] + x[9] + x[10] <= 1
 clq2: x[2] + x[3] + x[4] <= 1
 clq3: x[5] + x[6] + x[7] <= 1
 clq4: x[9] + x[7] + x[8] <= 1
 clq5: x[1] + x[2] <= 1
 clq6: x[3] + x[5] <= 1
 clq7: x[10] + x[8] <= 1
 clq8: x[4] + x[6] <= 1
 clq9: x[1] + x[3] <= 1
 clq10: x[2] + x[10] <= 1
 clq11: x[4] + x[5] <= 1
 clq12: x[6] + x[8] <= 1
Bounds
Binaries
 x[1] x[2] x[9] x[3] x[10] x[4] x[5] x[6] x[7] x[8]
End
