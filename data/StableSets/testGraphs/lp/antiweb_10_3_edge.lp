\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[2] + x[9] + x[3] + x[10] + x[4] + x[5] + x[6] + x[7] + x[8]
Subject To
 edge1: x[1] + x[2] <= 1
 edge2: x[1] + x[9] <= 1
 edge3: x[1] + x[3] <= 1
 edge4: x[1] + x[10] <= 1
 edge5: x[2] + x[4] <= 1
 edge6: x[2] + x[3] <= 1
 edge7: x[2] + x[10] <= 1
 edge8: x[9] + x[7] <= 1
 edge9: x[9] + x[8] <= 1
 edge10: x[9] + x[10] <= 1
 edge11: x[3] + x[4] <= 1
 edge12: x[3] + x[5] <= 1
 edge13: x[10] + x[8] <= 1
 edge14: x[4] + x[5] <= 1
 edge15: x[4] + x[6] <= 1
 edge16: x[5] + x[6] <= 1
 edge17: x[5] + x[7] <= 1
 edge18: x[6] + x[7] <= 1
 edge19: x[6] + x[8] <= 1
 edge20: x[7] + x[8] <= 1
Bounds
Binaries
 x[1] x[2] x[9] x[3] x[10] x[4] x[5] x[6] x[7] x[8]
End
