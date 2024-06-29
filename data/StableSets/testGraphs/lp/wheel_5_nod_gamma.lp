\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[3] + x[4] + x[5] + x[6]
Subject To
 nod1: 5 x[1] + x[2] + x[3] + x[4] + x[5] + x[6] <= 5
 nod2: x[1] + 3 x[2] + x[3] + x[6] <= 3
 nod3: x[1] + x[2] + 3 x[3] + x[4] <= 3
 nod4: x[1] + x[3] + 3 x[4] + x[5] <= 3
 nod5: x[1] + x[4] + 3 x[5] + x[6] <= 3
 nod6: x[1] + x[2] + x[5] + 3 x[6] <= 3
Bounds
Binaries
 x[1] x[2] x[3] x[4] x[5] x[6]
End
