\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[7] + x[3] + x[8] + x[4] + x[5] + x[6]
Subject To
 nod1: 4 x[1] + x[2] + x[7] + x[3] + x[8] <= 4
 nod2: x[1] + 4 x[2] + x[3] + x[8] + x[4] <= 4
 nod7: x[1] + 4 x[7] + x[8] + x[5] + x[6] <= 4
 nod3: x[1] + x[2] + 4 x[3] + x[4] + x[5] <= 4
 nod8: x[1] + x[2] + x[7] + 4 x[8] + x[6] <= 4
 nod4: x[2] + x[3] + 4 x[4] + x[5] + x[6] <= 4
 nod5: x[7] + x[3] + x[4] + 4 x[5] + x[6] <= 4
 nod6: x[7] + x[8] + x[4] + x[5] + 4 x[6] <= 4
Bounds
Binaries
 x[1] x[2] x[7] x[3] x[8] x[4] x[5] x[6]
End
