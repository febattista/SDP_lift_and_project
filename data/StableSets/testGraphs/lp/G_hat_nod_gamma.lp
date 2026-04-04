\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x8a0cf6aa34878dfe
Maximize
obj:   x[1] + x[2] + x[4] + x[5] + x[6] + x[3] + x[7] + x[8]
Subject To
 nod1: 4 x[1] + x[2] + x[4] + x[5] + x[6] <= 4
 nod2: x[1] + 4 x[2] + x[6] + x[3] + x[7] <= 4
 nod4: x[1] + 4 x[4] + x[5] + x[3] + x[8] <= 4
 nod5: x[1] + x[4] + 5 x[5] + x[6] + x[7] + x[8] <= 5
 nod6: x[1] + x[2] + x[5] + 5 x[6] + x[7] + x[8] <= 5
 nod3: x[2] + x[4] + 4 x[3] + x[7] + x[8] <= 4
 nod7: x[2] + x[5] + x[6] + x[3] + 5 x[7] + x[8] <= 5
 nod8: x[4] + x[5] + x[6] + x[3] + x[7] + 5 x[8] <= 5
Bounds
Binaries
 x[1] x[2] x[4] x[5] x[6] x[3] x[7] x[8]
End
