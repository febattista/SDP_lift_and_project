\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[7] + x[3] + x[8] + x[4] + x[5] + x[6]
Subject To
 nod1: 2 x[1] + x[2] + x[7] + x[3] + x[8] <= 2
 nod2: x[1] + 2 x[2] + x[3] + x[8] + x[4] <= 2
 nod7: x[1] + 2 x[7] + x[8] + x[5] + x[6] <= 2
 nod3: x[1] + x[2] + 2 x[3] + x[4] + x[5] <= 2
 nod8: x[1] + x[2] + x[7] + 2 x[8] + x[6] <= 2
 nod4: x[2] + x[3] + 2 x[4] + x[5] + x[6] <= 2
 nod5: x[7] + x[3] + x[4] + 2 x[5] + x[6] <= 2
 nod6: x[7] + x[8] + x[4] + x[5] + 2 x[6] <= 2
Bounds
Binaries
 x[1] x[2] x[7] x[3] x[8] x[4] x[5] x[6]
End
