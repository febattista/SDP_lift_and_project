\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[9] + x[3] + x[10] + x[4] + x[5] + x[6] + x[7] + x[8]
Subject To
 nod1: 2 x[1] + x[2] + x[9] + x[3] + x[10] <= 2
 nod2: x[1] + 2 x[2] + x[3] + x[10] + x[4] <= 2
 nod9: x[1] + 2 x[9] + x[10] + x[7] + x[8] <= 2
 nod3: x[1] + x[2] + 2 x[3] + x[4] + x[5] <= 2
 nod10: x[1] + x[2] + x[9] + 2 x[10] + x[8] <= 2
 nod4: x[2] + x[3] + 2 x[4] + x[5] + x[6] <= 2
 nod5: x[3] + x[4] + 2 x[5] + x[6] + x[7] <= 2
 nod6: x[4] + x[5] + 2 x[6] + x[7] + x[8] <= 2
 nod7: x[9] + x[5] + x[6] + 2 x[7] + x[8] <= 2
 nod8: x[9] + x[10] + x[6] + x[7] + 2 x[8] <= 2
Bounds
Binaries
 x[1] x[2] x[9] x[3] x[10] x[4] x[5] x[6] x[7] x[8]
End
