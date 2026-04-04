\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0xdd8608e257b6b071
Maximize
obj:   x[1] + x[2] + x[4] + x[5] + x[6] + x[3] + x[7] + x[8]
Subject To
 clq1: x[1] + x[4] + x[5] <= 1
 clq2: x[2] + x[3] + x[7] <= 1
 clq3: x[1] + x[2] + x[6] <= 1
 clq4: x[4] + x[3] + x[8] <= 1
 clq5: x[5] + x[6] + x[7] + x[8] <= 1
Bounds
Binaries
 x[1] x[2] x[4] x[5] x[6] x[3] x[7] x[8]
End
