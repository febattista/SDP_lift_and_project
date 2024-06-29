\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[9] + x[3] + x[10] + x[4] + x[5] + x[6] + x[7] + x[8]
Subject To
 nod1: 4 x[1] + x[2] + x[9] + x[3] + x[10] <= 4
 nod2: x[1] + 4 x[2] + x[3] + x[10] + x[4] <= 4
 nod9: x[1] + 4 x[9] + x[10] + x[7] + x[8] <= 4
 nod3: x[1] + x[2] + 4 x[3] + x[4] + x[5] <= 4
 nod10: x[1] + x[2] + x[9] + 4 x[10] + x[8] <= 4
 nod4: x[2] + x[3] + 4 x[4] + x[5] + x[6] <= 4
 nod5: x[3] + x[4] + 4 x[5] + x[6] + x[7] <= 4
 nod6: x[4] + x[5] + 4 x[6] + x[7] + x[8] <= 4
 nod7: x[9] + x[5] + x[6] + 4 x[7] + x[8] <= 4
 nod8: x[9] + x[10] + x[6] + x[7] + 4 x[8] <= 4
Bounds
Binaries
 x[1] x[2] x[9] x[3] x[10] x[4] x[5] x[6] x[7] x[8]
End
