\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[128] + x[129] + x[67] + x[70] + x[7] + x[8] + x[39] + x[89]
   + x[77] + x[143] + x[27] + x[115] + x[118] + x[23] + x[121] + x[87]
   + x[79] + x[94] + x[2] + x[32] + x[65] + x[135] + x[131] + x[155]
   + x[134] + x[72] + x[73] + x[139] + x[76] + x[13] + x[136] + x[148]
   + x[169] + x[111] + x[28] + x[93] + x[62] + x[3] + x[11] + x[140]
   + x[172] + x[153] + x[4] + x[130] + x[16] + x[24] + x[34] + x[37]
   + x[171] + x[175] + x[57] + x[58] + x[78] + x[88] + x[92] + x[97]
   + x[100] + x[103] + x[109] + x[125] + x[126] + x[127] + x[5] + x[160]
   + x[48] + x[166] + x[41] + x[108] + x[141] + x[44] + x[81] + x[146]
   + x[22] + x[119] + x[25] + x[29] + x[158] + x[6] + x[98] + x[71]
   + x[104] + x[42] + x[12] + x[80] + x[83] + x[117] + x[74] + x[96]
   + x[107] + x[75] + x[113] + x[147] + x[53] + x[56] + x[152] + x[91]
   + x[95] + x[133] + x[161] + x[50] + x[9] + x[123] + x[165] + x[35]
   + x[68] + x[38] + x[106] + x[112] + x[84] + x[120] + x[54] + x[60]
   + x[10] + x[114] + x[36] + x[26] + x[66] + x[116] + x[52] + x[20]
   + x[33] + x[149] + x[154] + x[30] + x[110] + x[124] + x[14] + x[64]
   + x[163] + x[18] + x[21] + x[122] + x[15] + x[90] + x[170] + x[19]
   + x[40] + x[82] + x[105] + x[63] + x[17] + x[102] + x[45] + x[46]
   + x[59] + x[132] + x[174] + x[138] + x[101] + x[151] + x[173] + x[99]
   + x[85] + x[86] + x[150] + x[49] + x[142] + x[61] + x[144] + x[145]
   + x[43] + x[55] + x[51] + x[164] + x[156] + x[137] + x[31] + x[69]
   + x[162] + x[159] + x[168] + x[167] + x[47] + x[157]
Subject To
 nod1: 11 x[1] + x[128] + x[129] + x[67] + x[70] + x[7] + x[8] + x[39]
   + x[89] + x[77] + x[143] + x[27] + x[115] + x[118] + x[23] + x[121]
   + x[87] + x[79] + x[94] <= 11
 nod128: x[1] + 12 x[128] + x[76] + x[153] + x[6] + x[80] + x[113] + x[95]
   + x[161] + x[50] + x[149] + x[30] + x[110] + x[64] + x[46] + x[132]
   + x[138] + x[51] <= 12
 nod129: x[1] + 14 x[129] + x[70] + x[115] + x[4] + x[16] + x[171] + x[97]
   + x[100] + x[160] + x[108] + x[81] + x[119] + x[158] + x[104] + x[80]
   + x[91] + x[116] + x[14] + x[15] + x[19] + x[144] + x[55] <= 14
 nod67: x[1] + 15 x[67] + x[115] + x[118] + x[94] + x[76] + x[28] + x[175]
   + x[5] + x[147] + x[56] + x[112] + x[20] + x[149] + x[154] + x[110]
   + x[122] + x[15] + x[86] + x[150] + x[145] + x[55] + x[31] <= 15
 nod70: x[1] + x[129] + 11 x[70] + x[8] + x[118] + x[134] + x[139] + x[169]
   + x[5] + x[160] + x[113] + x[147] + x[9] + x[106] + x[114] + x[64]
   + x[170] + x[102] + x[43] + x[168] + x[47] <= 11
 nod7: x[1] + 11 x[7] + x[134] + x[76] + x[78] + x[88] + x[109] + x[25]
   + x[29] + x[96] + x[107] + x[75] + x[113] + x[147] + x[53] + x[56]
   + x[152] + x[91] + x[95] <= 11
 nod8: x[1] + x[70] + 7 x[8] + x[39] + x[153] + x[100] + x[81] + x[152]
   + x[133] + x[161] + x[50] + x[9] + x[123] + x[165] <= 7
 nod39: x[1] + x[8] + 13 x[39] + x[2] + x[65] + x[136] + x[62] + x[11]
   + x[16] + x[37] + x[171] + x[160] + x[81] + x[71] + x[74] + x[161]
   + x[123] + x[165] + x[116] + x[20] + x[30] + x[150] + x[49] + x[144]
   + x[167] + x[47] <= 13
 nod89: x[1] + 14 x[89] + x[72] + x[148] + x[172] + x[37] + x[126] + x[119]
   + x[71] + x[113] + x[152] + x[106] + x[116] + x[52] + x[154] + x[18]
   + x[82] + x[150] + x[31] + x[168] <= 14
 nod77: x[1] + 10 x[77] + x[27] + x[140] + x[92] + x[103] + x[125] + x[81]
   + x[74] + x[95] + x[54] + x[82] + x[173] + x[145] + x[156] + x[162]
   <= 10
 nod143: x[1] + 9 x[143] + x[27] + x[23] + x[87] + x[4] + x[37] + x[127]
   + x[166] + x[98] + x[83] + x[74] + x[110] + x[14] + x[90] + x[101]
   + x[157] <= 9
 nod27: x[1] + x[77] + x[143] + 9 x[27] + x[172] + x[37] + x[92] + x[125]
   + x[141] + x[74] + x[53] + x[84] + x[66] + x[82] <= 9
 nod115: x[1] + x[129] + x[67] + 10 x[115] + x[13] + x[34] + x[125] + x[41]
   + x[107] + x[161] + x[116] + x[149] + x[154] + x[30] + x[122] + x[90]
   + x[132] + x[85] <= 10
 nod118: x[1] + x[67] + x[70] + 12 x[118] + x[79] + x[135] + x[139]
   + x[140] + x[146] + x[158] + x[71] + x[83] + x[107] + x[152] + x[84]
   + x[10] + x[105] + x[102] + x[137] <= 12
 nod23: x[1] + x[143] + 14 x[23] + x[87] + x[109] + x[22] + x[104] + x[96]
   + x[152] + x[9] + x[165] + x[68] + x[84] + x[54] + x[116] + x[52]
   + x[174] + x[151] + x[86] + x[150] + x[49] + x[144] + x[145] + x[43]
   + x[55] <= 14
 nod121: x[1] + 11 x[121] + x[32] + x[134] + x[28] + x[93] + x[4] + x[109]
   + x[29] + x[152] + x[161] + x[10] + x[20] + x[90] + x[174] + x[145]
   <= 11
 nod87: x[1] + x[143] + x[23] + 10 x[87] + x[155] + x[92] + x[127] + x[166]
   + x[41] + x[81] + x[83] + x[113] + x[35] + x[21] + x[82] + x[85] <= 10
 nod79: x[1] + x[118] + 14 x[79] + x[2] + x[32] + x[3] + x[11] + x[37]
   + x[97] + x[100] + x[126] + x[104] + x[42] + x[152] + x[95] + x[9]
   + x[165] + x[36] + x[33] + x[14] + x[122] + x[15] + x[170] + x[40]
   + x[138] + x[51] <= 14
 nod94: x[1] + x[67] + 13 x[94] + x[139] + x[3] + x[48] + x[166] + x[141]
   + x[44] + x[71] + x[83] + x[10] + x[52] + x[21] + x[170] + x[46] + x[59]
   + x[150] + x[49] + x[43] + x[157] <= 13
 nod2: x[39] + x[79] + 13 x[2] + x[32] + x[65] + x[135] + x[131] + x[155]
   + x[134] + x[72] + x[73] + x[139] + x[76] + x[13] + x[136] + x[148]
   + x[169] + x[111] + x[28] + x[93] + x[62] <= 13
 nod32: x[121] + x[79] + x[2] + 8 x[32] + x[72] + x[76] + x[62] + x[3]
   + x[37] + x[152] + x[21] + x[63] + x[46] + x[173] + x[162] + x[159]
   <= 8
 nod65: x[39] + x[2] + 11 x[65] + x[141] + x[146] + x[22] + x[75] + x[95]
   + x[124] + x[40] + x[82] + x[46] + x[59] + x[101] + x[61] + x[43]
   + x[156] <= 11
 nod135: x[118] + x[2] + 12 x[135] + x[171] + x[88] + x[22] + x[42] + x[91]
   + x[133] + x[165] + x[154] + x[30] + x[86] + x[61] + x[156] + x[31]
   + x[47] <= 12
 nod131: x[2] + 12 x[131] + x[111] + x[3] + x[125] + x[160] + x[119]
   + x[25] + x[29] + x[98] + x[117] + x[152] + x[68] + x[124] + x[163]
   + x[18] + x[132] + x[51] + x[137] <= 12
 nod155: x[87] + x[2] + 12 x[155] + x[148] + x[92] + x[5] + x[81] + x[113]
   + x[133] + x[123] + x[106] + x[112] + x[26] + x[33] + x[30] + x[19]
   + x[82] + x[132] + x[101] + x[151] + x[156] + x[159] <= 12
 nod134: x[70] + x[7] + x[121] + x[2] + 11 x[134] + x[111] + x[97] + x[12]
   + x[83] + x[75] + x[165] + x[33] + x[110] + x[21] + x[15] + x[63]
   + x[132] + x[145] <= 11
 nod72: x[89] + x[2] + x[32] + 14 x[72] + x[13] + x[28] + x[62] + x[3]
   + x[171] + x[175] + x[58] + x[100] + x[166] + x[141] + x[81] + x[120]
   + x[14] + x[63] + x[61] + x[145] + x[51] + x[162] <= 14
 nod73: x[2] + 12 x[73] + x[169] + x[28] + x[3] + x[175] + x[103] + x[161]
   + x[9] + x[110] + x[18] + x[105] + x[46] + x[174] + x[101] + x[99]
   + x[150] + x[142] + x[162] <= 12
 nod139: x[70] + x[118] + x[94] + x[2] + 11 x[139] + x[4] + x[130] + x[34]
   + x[22] + x[74] + x[147] + x[95] + x[123] + x[165] + x[21] + x[15]
   + x[90] + x[82] + x[174] <= 11
 nod76: x[128] + x[67] + x[7] + x[2] + x[32] + 15 x[76] + x[100] + x[103]
   + x[48] + x[108] + x[25] + x[12] + x[117] + x[50] + x[35] + x[112]
   + x[84] + x[60] + x[114] + x[52] + x[18] + x[122] + x[45] + x[174]
   + x[145] + x[167] <= 15
 nod13: x[115] + x[2] + x[72] + 8 x[13] + x[57] + x[125] + x[83] + x[107]
   + x[75] + x[113] + x[35] + x[52] + x[110] + x[124] <= 8
 nod136: x[39] + x[2] + 14 x[136] + x[28] + x[130] + x[34] + x[100]
   + x[103] + x[126] + x[160] + x[48] + x[166] + x[22] + x[104] + x[12]
   + x[96] + x[91] + x[165] + x[54] + x[36] + x[33] + x[105] + x[17]
   + x[151] + x[31] <= 14
 nod148: x[89] + x[2] + x[155] + 8 x[148] + x[125] + x[160] + x[6] + x[113]
   + x[68] + x[26] + x[116] + x[122] + x[43] <= 8
 nod169: x[70] + x[2] + x[73] + 11 x[169] + x[88] + x[152] + x[35] + x[84]
   + x[36] + x[15] + x[170] + x[63] + x[45] + x[138] + x[101] + x[49]
   <= 11
 nod111: x[2] + x[131] + x[134] + 11 x[111] + x[24] + x[41] + x[29] + x[12]
   + x[147] + x[91] + x[165] + x[10] + x[114] + x[116] + x[30] + x[15]
   + x[45] + x[49] + x[144] + x[157] <= 11
 nod28: x[67] + x[121] + x[2] + x[72] + x[73] + x[136] + 11 x[28] + x[34]
   + x[171] + x[58] + x[97] + x[160] + x[158] + x[12] + x[75] + x[165]
   + x[66] + x[17] + x[86] <= 11
 nod93: x[121] + x[2] + 12 x[93] + x[153] + x[4] + x[24] + x[34] + x[57]
   + x[78] + x[108] + x[83] + x[74] + x[147] + x[9] + x[114] + x[110]
   + x[122] + x[47] <= 12
 nod62: x[39] + x[2] + x[32] + x[72] + 10 x[62] + x[166] + x[81] + x[25]
   + x[71] + x[107] + x[161] + x[60] + x[116] + x[122] + x[82] + x[63]
   + x[102] + x[45] + x[150] + x[49] + x[157] <= 10
 nod3: x[79] + x[94] + x[32] + x[131] + x[72] + x[73] + 8 x[3] + x[11]
   + x[140] + x[172] + x[153] <= 8
 nod11: x[39] + x[79] + x[3] + 9 x[11] + x[172] + x[57] + x[58] + x[29]
   + x[98] + x[42] + x[91] + x[60] + x[66] + x[116] + x[52] + x[20] <= 9
 nod140: x[77] + x[118] + x[3] + 8 x[140] + x[103] + x[125] + x[83]
   + x[107] + x[50] + x[68] + x[112] + x[40] + x[162] <= 8
 nod172: x[89] + x[27] + x[3] + x[11] + 10 x[172] + x[108] + x[119] + x[71]
   + x[91] + x[9] + x[122] + x[45] + x[51] + x[168] <= 10
 nod153: x[128] + x[8] + x[93] + x[3] + 12 x[153] + x[5] + x[166] + x[6]
   + x[107] + x[152] + x[84] + x[60] + x[114] + x[14] + x[63] + x[174]
   + x[85] <= 12
 nod4: x[129] + x[143] + x[121] + x[139] + x[93] + 14 x[4] + x[130] + x[16]
   + x[24] + x[34] + x[37] + x[171] + x[175] + x[57] + x[58] + x[78]
   + x[88] + x[92] + x[97] + x[100] + x[103] + x[109] + x[125] + x[126]
   + x[127] <= 14
 nod130: x[139] + x[136] + x[4] + 11 x[130] + x[103] + x[80] + x[74]
   + x[113] + x[152] + x[120] + x[10] + x[66] + x[163] + x[99] + x[49]
   + x[159] <= 11
 nod16: x[129] + x[39] + x[4] + 9 x[16] + x[97] + x[5] + x[98] + x[147]
   + x[133] + x[9] + x[60] + x[36] + x[124] + x[40] + x[82] + x[105]
   + x[63] <= 9
 nod24: x[111] + x[93] + x[4] + 10 x[24] + x[125] + x[160] + x[81] + x[91]
   + x[112] + x[124] + x[40] + x[82] + x[59] + x[99] + x[51] + x[164]
   <= 10
 nod34: x[115] + x[139] + x[136] + x[28] + x[93] + x[4] + 11 x[34] + x[96]
   + x[53] + x[91] + x[20] + x[33] + x[19] + x[17] + x[55] + x[51] + x[69]
   <= 11
 nod37: x[39] + x[89] + x[143] + x[27] + x[79] + x[32] + x[4] + 11 x[37]
   + x[126] + x[160] + x[107] + x[50] + x[123] + x[68] + x[112] + x[120]
   + x[114] + x[116] + x[52] + x[86] + x[167] <= 11
 nod171: x[129] + x[39] + x[135] + x[72] + x[28] + x[4] + 5 x[171] + x[17]
   + x[102] + x[49] <= 5
 nod175: x[67] + x[72] + x[73] + x[4] + 12 x[175] + x[22] + x[119] + x[75]
   + x[60] + x[36] + x[52] + x[149] + x[163] + x[45] + x[69] <= 12
 nod57: x[13] + x[93] + x[11] + x[4] + 10 x[57] + x[103] + x[125] + x[48]
   + x[108] + x[141] + x[50] + x[36] + x[110] + x[19] + x[151] + x[142]
   + x[43] <= 10
 nod58: x[72] + x[28] + x[11] + x[4] + 10 x[58] + x[78] + x[126] + x[165]
   + x[114] + x[52] + x[154] + x[110] + x[163] + x[15] + x[132] + x[137]
   + x[31] + x[167] <= 10
 nod78: x[7] + x[93] + x[4] + x[58] + 10 x[78] + x[109] + x[81] + x[74]
   + x[147] + x[95] + x[36] + x[66] + x[30] + x[14] + x[151] + x[49]
   + x[31] <= 10
 nod88: x[7] + x[135] + x[169] + x[4] + 10 x[88] + x[109] + x[81] + x[158]
   + x[6] + x[123] + x[165] + x[124] + x[150] + x[55] <= 10
 nod92: x[77] + x[27] + x[87] + x[155] + x[4] + 8 x[92] + x[103] + x[125]
   + x[41] + x[12] + x[165] + x[163] + x[101] + x[86] + x[137] + x[167]
   <= 8
 nod97: x[129] + x[79] + x[134] + x[28] + x[4] + x[16] + 8 x[97] + x[160]
   + x[12] + x[74] + x[53] + x[154] + x[163] + x[151] + x[51] <= 8
 nod100: x[129] + x[8] + x[79] + x[72] + x[76] + x[136] + x[4] + 14 x[100]
   + x[125] + x[126] + x[81] + x[12] + x[152] + x[95] + x[112] + x[60]
   + x[10] + x[116] + x[14] + x[122] + x[90] + x[19] + x[45] + x[46]
   + x[132] + x[55] + x[137] <= 14
 nod103: x[77] + x[73] + x[76] + x[136] + x[140] + x[4] + x[130] + x[57]
   + x[92] + 14 x[103] + x[5] + x[25] + x[42] + x[20] + x[154] + x[30]
   + x[64] + x[21] + x[90] + x[101] + x[145] + x[69] + x[162] <= 14
 nod109: x[7] + x[23] + x[121] + x[4] + x[78] + x[88] + 9 x[109] + x[126]
   + x[127] + x[41] + x[29] + x[71] + x[117] + x[123] + x[163] + x[122]
   + x[144] + x[145] <= 9
 nod125: x[77] + x[27] + x[115] + x[131] + x[13] + x[148] + x[140] + x[4]
   + x[24] + x[57] + x[92] + x[100] + 15 x[125] + x[104] + x[107] + x[161]
   + x[123] + x[33] + x[90] + x[82] + x[17] + x[144] + x[145] + x[55]
   + x[157] <= 15
 nod126: x[89] + x[79] + x[136] + x[4] + x[37] + x[58] + x[100] + x[109]
   + 11 x[126] + x[44] + x[74] + x[96] + x[75] + x[161] + x[165] + x[112]
   + x[144] + x[137] <= 11
 nod127: x[143] + x[87] + x[4] + x[109] + 7 x[127] + x[71] + x[133]
   + x[161] + x[50] + x[112] + x[15] + x[162] <= 7
 nod5: x[67] + x[70] + x[155] + x[153] + x[16] + x[103] + 14 x[5] + x[160]
   + x[48] + x[166] + x[41] + x[108] + x[141] + x[44] + x[81] + x[146]
   + x[22] + x[119] + x[25] + x[29] + x[158] <= 14
 nod160: x[129] + x[70] + x[39] + x[131] + x[136] + x[148] + x[28] + x[24]
   + x[37] + x[97] + x[5] + 13 x[160] + x[48] + x[6] + x[74] + x[161]
   + x[36] + x[122] + x[170] + x[105] + x[63] + x[101] + x[43] + x[55]
   + x[51] + x[164] + x[159] <= 13
 nod48: x[94] + x[76] + x[136] + x[57] + x[5] + x[160] + 8 x[48] + x[25]
   + x[75] + x[120] + x[18] + x[19] + x[31] <= 8
 nod166: x[143] + x[87] + x[94] + x[72] + x[136] + x[62] + x[153] + x[5]
   + 14 x[166] + x[44] + x[25] + x[29] + x[95] + x[9] + x[54] + x[66]
   + x[116] + x[20] + x[122] + x[102] + x[43] + x[167] <= 14
 nod41: x[115] + x[87] + x[111] + x[92] + x[109] + x[5] + 12 x[41] + x[71]
   + x[107] + x[112] + x[60] + x[26] + x[52] + x[20] + x[163] + x[122]
   + x[170] + x[151] + x[85] + x[49] + x[51] + x[157] <= 12
 nod108: x[129] + x[76] + x[93] + x[172] + x[57] + x[5] + 11 x[108] + x[81]
   + x[113] + x[50] + x[149] + x[170] + x[40] + x[105] + x[46] + x[51]
   + x[157] <= 11
 nod141: x[27] + x[94] + x[65] + x[72] + x[57] + x[5] + 11 x[141] + x[44]
   + x[25] + x[104] + x[36] + x[14] + x[122] + x[15] + x[90] + x[19]
   + x[59] + x[132] + x[174] <= 11
 nod44: x[94] + x[126] + x[5] + x[166] + x[141] + 12 x[44] + x[158] + x[71]
   + x[12] + x[83] + x[96] + x[152] + x[133] + x[165] + x[114] + x[40]
   + x[59] + x[174] + x[101] + x[173] + x[86] + x[51] + x[164] + x[157]
   <= 12
 nod81: x[129] + x[8] + x[39] + x[77] + x[87] + x[155] + x[72] + x[62]
   + x[24] + x[78] + x[88] + x[100] + x[5] + x[108] + 14 x[81] + x[22]
   + x[98] + x[147] + x[56] + x[38] + x[14] + x[40] + x[82] + x[105]
   + x[138] + x[150] + x[137] <= 14
 nod146: x[118] + x[65] + x[5] + 8 x[146] + x[29] + x[152] + x[35] + x[30]
   + x[15] + x[82] + x[59] + x[51] + x[69] <= 8
 nod22: x[23] + x[65] + x[135] + x[139] + x[136] + x[175] + x[5] + x[81]
   + 11 x[22] + x[42] + x[80] + x[9] + x[120] + x[33] + x[149] + x[102]
   + x[174] + x[85] + x[49] + x[142] + x[61] <= 11
 nod119: x[129] + x[89] + x[131] + x[172] + x[175] + x[5] + 12 x[119]
   + x[161] + x[38] + x[84] + x[54] + x[66] + x[163] + x[132] + x[138]
   + x[101] <= 12
 nod25: x[7] + x[131] + x[76] + x[62] + x[103] + x[5] + x[48] + x[166]
   + x[141] + 12 x[25] + x[80] + x[83] + x[112] + x[54] + x[114] + x[149]
   + x[30] + x[64] + x[55] + x[164] + x[156] <= 12
 nod29: x[7] + x[121] + x[131] + x[111] + x[11] + x[109] + x[5] + x[166]
   + x[146] + 13 x[29] + x[75] + x[53] + x[95] + x[133] + x[38] + x[120]
   + x[19] + x[17] + x[85] + x[86] + x[61] + x[156] <= 13
 nod158: x[129] + x[118] + x[28] + x[88] + x[5] + x[44] + 12 x[158] + x[71]
   + x[107] + x[113] + x[95] + x[10] + x[149] + x[110] + x[40] + x[82]
   + x[47] + x[157] <= 12
 nod6: x[128] + x[148] + x[153] + x[88] + x[160] + 9 x[6] + x[98] + x[71]
   + x[104] + x[42] + x[12] + x[80] + x[83] + x[117] + x[74] <= 9
 nod98: x[143] + x[131] + x[11] + x[16] + x[81] + x[6] + 8 x[98] + x[42]
   + x[113] + x[56] + x[133] + x[66] + x[116] + x[124] <= 8
 nod71: x[39] + x[89] + x[118] + x[94] + x[62] + x[172] + x[109] + x[127]
   + x[41] + x[44] + x[158] + x[6] + 9 x[71] + x[42] + x[116] + x[20]
   + x[15] + x[151] + x[150] + x[144] <= 9
 nod104: x[129] + x[23] + x[79] + x[136] + x[125] + x[141] + x[6]
   + 11 x[104] + x[152] + x[66] + x[124] + x[18] + x[122] + x[63] + x[59]
   + x[174] + x[156] <= 11
 nod42: x[79] + x[135] + x[11] + x[103] + x[22] + x[6] + x[98] + x[71]
   + 10 x[42] + x[12] + x[152] + x[36] + x[149] + x[19] + x[63] + x[174]
   + x[85] + x[144] + x[31] <= 10
 nod12: x[134] + x[76] + x[136] + x[111] + x[28] + x[92] + x[97] + x[100]
   + x[44] + x[6] + x[42] + 9 x[12] + x[83] + x[152] + x[112] + x[60]
   + x[33] + x[149] + x[154] + x[30] <= 9
 nod80: x[128] + x[129] + x[130] + x[22] + x[25] + x[6] + 11 x[80] + x[74]
   + x[10] + x[40] + x[101] + x[151] + x[145] + x[159] <= 11
 nod83: x[143] + x[118] + x[87] + x[94] + x[134] + x[13] + x[93] + x[140]
   + x[44] + x[25] + x[6] + x[12] + 14 x[83] + x[9] + x[68] + x[54]
   + x[110] + x[122] + x[170] + x[102] + x[138] + x[144] + x[55] + x[69]
   <= 14
 nod117: x[131] + x[76] + x[109] + x[6] + 7 x[117] + x[113] + x[165]
   + x[52] + x[59] + x[43] + x[137] + x[162] <= 7
 nod74: x[39] + x[77] + x[143] + x[27] + x[139] + x[93] + x[130] + x[78]
   + x[97] + x[126] + x[160] + x[6] + x[80] + 12 x[74] + x[53] + x[123]
   + x[60] + x[105] + x[102] + x[59] <= 12
 nod96: x[7] + x[23] + x[136] + x[34] + x[126] + x[44] + 9 x[96] + x[60]
   + x[114] + x[26] + x[163] + x[174] + x[101] + x[49] + x[157] <= 9
 nod107: x[7] + x[115] + x[118] + x[13] + x[62] + x[140] + x[153] + x[37]
   + x[125] + x[41] + x[158] + 12 x[107] + x[56] + x[120] + x[54] + x[114]
   + x[124] + x[150] + x[61] + x[145] + x[164] + x[168] <= 12
 nod75: x[7] + x[65] + x[134] + x[13] + x[28] + x[175] + x[126] + x[48]
   + x[29] + 13 x[75] + x[33] + x[64] + x[21] + x[122] + x[99] + x[164]
   + x[156] <= 13
 nod113: x[128] + x[70] + x[7] + x[89] + x[87] + x[155] + x[13] + x[148]
   + x[130] + x[108] + x[158] + x[98] + x[117] + 15 x[113] + x[161]
   + x[165] + x[10] + x[26] + x[20] + x[64] + x[19] + x[138] + x[43]
   + x[55] + x[47] <= 15
 nod147: x[67] + x[70] + x[7] + x[139] + x[111] + x[93] + x[16] + x[78]
   + x[81] + 9 x[147] + x[152] + x[82] + x[173] + x[150] + x[144] + x[164]
   <= 9
 nod53: x[7] + x[27] + x[34] + x[97] + x[29] + x[74] + 9 x[53] + x[84]
   + x[36] + x[85] + x[49] + x[61] + x[144] + x[145] <= 9
 nod56: x[67] + x[7] + x[81] + x[98] + x[107] + 6 x[56] + x[9] + x[68]
   + x[60] + x[52] + x[45] <= 6
 nod152: x[7] + x[8] + x[89] + x[118] + x[23] + x[121] + x[79] + x[32]
   + x[131] + x[169] + x[153] + x[130] + x[100] + x[44] + x[146] + x[104]
   + x[42] + x[12] + x[147] + 13 x[152] + x[95] + x[106] + x[33] + x[149]
   + x[170] + x[40] + x[63] + x[51] <= 13
 nod91: x[129] + x[7] + x[135] + x[136] + x[111] + x[11] + x[172] + x[24]
   + x[34] + 13 x[91] + x[114] + x[66] + x[52] + x[20] + x[163] + x[174]
   + x[31] + x[69] <= 13
 nod95: x[128] + x[7] + x[77] + x[79] + x[65] + x[139] + x[78] + x[100]
   + x[166] + x[29] + x[158] + x[152] + 12 x[95] + x[18] + x[59] + x[174]
   + x[49] + x[142] + x[144] <= 12
 nod133: x[8] + x[135] + x[155] + x[16] + x[127] + x[44] + x[29] + x[98]
   + 12 x[133] + x[50] + x[84] + x[116] + x[14] + x[142] + x[137] <= 12
 nod161: x[128] + x[8] + x[39] + x[115] + x[121] + x[73] + x[62] + x[125]
   + x[126] + x[127] + x[160] + x[119] + x[113] + 13 x[161] + x[35] + x[30]
   + x[145] + x[51] + x[159] + x[47] <= 13
 nod50: x[128] + x[8] + x[76] + x[140] + x[37] + x[57] + x[127] + x[108]
   + x[133] + 11 x[50] + x[120] + x[66] + x[99] + x[86] + x[31] + x[167]
   <= 11
 nod9: x[70] + x[8] + x[23] + x[79] + x[73] + x[93] + x[172] + x[16]
   + x[166] + x[22] + x[83] + x[56] + 14 x[9] + x[35] + x[68] + x[38]
   + x[106] + x[112] + x[84] + x[120] + x[54] + x[60] <= 14
 nod123: x[8] + x[39] + x[155] + x[139] + x[37] + x[88] + x[109] + x[125]
   + x[74] + 9 x[123] + x[165] + x[66] + x[46] + x[132] + x[173] <= 9
 nod165: x[8] + x[39] + x[23] + x[79] + x[135] + x[134] + x[139] + x[136]
   + x[111] + x[28] + x[58] + x[88] + x[92] + x[126] + x[44] + x[117]
   + x[113] + x[123] + 13 x[165] + x[84] + x[149] + x[59] + x[86] + x[144]
   + x[55] <= 13
 nod35: x[87] + x[76] + x[13] + x[169] + x[146] + x[161] + x[9] + 11 x[35]
   + x[66] + x[52] + x[163] + x[99] + x[61] <= 11
 nod68: x[23] + x[131] + x[148] + x[140] + x[37] + x[83] + x[56] + x[9]
   + 11 x[68] + x[60] + x[36] + x[90] + x[138] + x[151] + x[86] + x[49]
   + x[164] + x[159] <= 11
 nod38: x[81] + x[119] + x[29] + x[9] + 6 x[38] + x[59] + x[85] + x[156]
   <= 6
 nod106: x[70] + x[89] + x[155] + x[152] + x[9] + 7 x[106] + x[154] + x[64]
   + x[18] + x[170] + x[40] + x[43] <= 7
 nod112: x[67] + x[155] + x[76] + x[140] + x[24] + x[37] + x[100] + x[126]
   + x[127] + x[41] + x[25] + x[12] + x[9] + 11 x[112] + x[54] + x[173]
   + x[150] <= 11
 nod84: x[27] + x[118] + x[23] + x[76] + x[169] + x[153] + x[119] + x[53]
   + x[133] + x[9] + x[165] + 10 x[84] + x[114] + x[149] + x[168] <= 10
 nod120: x[72] + x[130] + x[37] + x[48] + x[22] + x[29] + x[107] + x[50]
   + x[9] + 13 x[120] + x[20] + x[33] + x[149] + x[138] + x[151] + x[85]
   + x[43] + x[55] + x[156] + x[137] <= 13
 nod54: x[77] + x[23] + x[136] + x[166] + x[119] + x[25] + x[83] + x[107]
   + x[9] + x[112] + 11 x[54] + x[10] + x[114] + x[154] + x[122] + x[63]
   + x[49] + x[157] <= 11
 nod60: x[76] + x[62] + x[11] + x[153] + x[16] + x[175] + x[100] + x[41]
   + x[12] + x[74] + x[96] + x[56] + x[9] + x[68] + 13 x[60] + x[33]
   + x[110] + x[18] + x[150] + x[55] + x[157] <= 13
 nod10: x[118] + x[121] + x[94] + x[111] + x[130] + x[100] + x[158] + x[80]
   + x[113] + x[54] + 9 x[10] + x[114] + x[36] + x[26] <= 9
 nod114: x[70] + x[76] + x[111] + x[93] + x[153] + x[37] + x[58] + x[44]
   + x[25] + x[96] + x[107] + x[91] + x[84] + x[54] + x[10] + 10 x[114]
   + x[116] + x[154] + x[64] + x[31] + x[168] <= 10
 nod36: x[79] + x[136] + x[169] + x[16] + x[175] + x[57] + x[78] + x[160]
   + x[141] + x[42] + x[53] + x[68] + x[10] + 14 x[36] + x[116] + x[149]
   + x[124] + x[122] + x[138] + x[101] + x[99] + x[51] <= 14
 nod26: x[155] + x[148] + x[41] + x[96] + x[113] + x[10] + 12 x[26]
   + x[154] + x[105] + x[59] + x[138] + x[101] + x[142] + x[61] + x[55]
   + x[137] <= 12
 nod66: x[27] + x[28] + x[11] + x[130] + x[78] + x[166] + x[119] + x[98]
   + x[104] + x[91] + x[50] + x[123] + x[35] + 14 x[66] + x[102] + x[132]
   + x[145] + x[156] + x[157] <= 14
 nod116: x[129] + x[39] + x[89] + x[115] + x[23] + x[148] + x[111] + x[62]
   + x[11] + x[37] + x[100] + x[166] + x[98] + x[71] + x[133] + x[114]
   + x[36] + 14 x[116] + x[40] + x[63] + x[46] + x[101] + x[159] + x[157]
   <= 14
 nod52: x[89] + x[23] + x[94] + x[76] + x[13] + x[11] + x[37] + x[175]
   + x[58] + x[41] + x[117] + x[56] + x[91] + x[35] + 12 x[52] + x[21]
   + x[40] + x[102] + x[144] <= 12
 nod20: x[67] + x[39] + x[121] + x[11] + x[34] + x[103] + x[166] + x[41]
   + x[71] + x[113] + x[91] + x[120] + 12 x[20] + x[170] + x[151] + x[173]
   + x[99] + x[85] + x[86] <= 12
 nod33: x[79] + x[155] + x[134] + x[136] + x[34] + x[125] + x[22] + x[12]
   + x[75] + x[152] + x[120] + x[60] + 9 x[33] + x[149] + x[163] + x[40]
   + x[82] + x[105] + x[132] + x[168] <= 9
 nod149: x[128] + x[67] + x[115] + x[175] + x[108] + x[22] + x[25] + x[158]
   + x[42] + x[12] + x[152] + x[165] + x[84] + x[120] + x[36] + x[33]
   + 14 x[149] + x[154] + x[18] + x[122] + x[173] + x[144] + x[43] + x[51]
   + x[167] <= 14
 nod154: x[67] + x[89] + x[115] + x[135] + x[58] + x[97] + x[103] + x[12]
   + x[106] + x[54] + x[114] + x[26] + x[149] + 11 x[154] + x[19] + x[173]
   + x[43] <= 11
 nod30: x[128] + x[39] + x[115] + x[135] + x[155] + x[111] + x[78] + x[103]
   + x[146] + x[25] + x[12] + x[161] + 12 x[30] + x[64] + x[19] + x[40]
   + x[17] + x[174] <= 12
 nod110: x[128] + x[67] + x[143] + x[134] + x[73] + x[13] + x[93] + x[57]
   + x[58] + x[158] + x[83] + x[60] + 13 x[110] + x[15] + x[19] + x[59]
   + x[174] + x[142] + x[145] + x[162] <= 13
 nod124: x[65] + x[131] + x[13] + x[16] + x[24] + x[88] + x[98] + x[104]
   + x[107] + x[36] + 10 x[124] + x[90] + x[45] + x[132] <= 10
 nod14: x[129] + x[143] + x[79] + x[72] + x[153] + x[78] + x[100] + x[141]
   + x[81] + x[133] + 11 x[14] + x[64] + x[163] + x[18] + x[21] + x[122]
   + x[15] + x[90] <= 11
 nod64: x[128] + x[70] + x[103] + x[25] + x[75] + x[113] + x[106] + x[114]
   + x[30] + x[14] + 8 x[64] + x[82] + x[137] + x[157] <= 8
 nod163: x[131] + x[130] + x[175] + x[58] + x[92] + x[97] + x[109] + x[41]
   + x[119] + x[96] + x[91] + x[35] + x[33] + x[14] + 14 x[163] + x[105]
   + x[45] + x[85] + x[43] + x[31] <= 14
 nod18: x[89] + x[131] + x[73] + x[76] + x[48] + x[104] + x[95] + x[106]
   + x[60] + x[149] + x[14] + 10 x[18] + x[19] + x[132] + x[174] + x[138]
   <= 10
 nod21: x[87] + x[94] + x[32] + x[134] + x[139] + x[103] + x[75] + x[52]
   + x[14] + 11 x[21] + x[40] + x[105] + x[101] + x[151] + x[173] + x[86]
   + x[150] <= 11
 nod122: x[67] + x[115] + x[79] + x[76] + x[148] + x[93] + x[62] + x[172]
   + x[100] + x[109] + x[160] + x[166] + x[41] + x[141] + x[104] + x[83]
   + x[75] + x[54] + x[36] + x[149] + x[14] + 15 x[122] + x[19] + x[49]
   + x[142] + x[162] <= 15
 nod15: x[129] + x[67] + x[79] + x[134] + x[139] + x[169] + x[111] + x[58]
   + x[127] + x[141] + x[146] + x[71] + x[110] + x[14] + 10 x[15] + x[170]
   + x[19] <= 10
 nod90: x[143] + x[115] + x[121] + x[139] + x[100] + x[103] + x[125]
   + x[141] + x[68] + x[124] + x[14] + 10 x[90] + x[43] <= 10
 nod170: x[70] + x[79] + x[94] + x[169] + x[160] + x[41] + x[108] + x[83]
   + x[152] + x[106] + x[20] + x[15] + 7 x[170] + x[101] + x[51] <= 7
 nod19: x[129] + x[155] + x[34] + x[57] + x[100] + x[48] + x[141] + x[29]
   + x[42] + x[113] + x[154] + x[30] + x[110] + x[18] + x[122] + x[15]
   + 12 x[19] + x[101] + x[151] <= 12
 nod40: x[79] + x[65] + x[140] + x[16] + x[24] + x[108] + x[44] + x[81]
   + x[158] + x[80] + x[152] + x[106] + x[116] + x[52] + x[33] + x[30]
   + x[21] + 13 x[40] + x[105] + x[69] <= 13
 nod82: x[89] + x[77] + x[27] + x[87] + x[65] + x[155] + x[139] + x[62]
   + x[16] + x[24] + x[125] + x[81] + x[146] + x[158] + x[147] + x[33]
   + x[64] + 14 x[82] + x[63] + x[102] + x[45] + x[138] + x[51] + x[162]
   + x[47] + x[157] <= 14
 nod105: x[118] + x[73] + x[136] + x[16] + x[160] + x[108] + x[81] + x[74]
   + x[26] + x[33] + x[163] + x[21] + x[40] + 9 x[105] + x[46] <= 9
 nod63: x[32] + x[134] + x[72] + x[169] + x[62] + x[153] + x[16] + x[160]
   + x[104] + x[42] + x[152] + x[54] + x[116] + x[82] + 11 x[63] + x[45]
   + x[168] <= 11
 nod17: x[136] + x[28] + x[34] + x[171] + x[125] + x[29] + x[30] + 8 x[17]
   + x[102] + x[45] + x[46] + x[59] <= 8
 nod102: x[70] + x[118] + x[62] + x[171] + x[166] + x[22] + x[83] + x[74]
   + x[66] + x[52] + x[82] + x[17] + 9 x[102] + x[150] + x[61] <= 9
 nod45: x[76] + x[169] + x[111] + x[62] + x[172] + x[175] + x[100] + x[56]
   + x[124] + x[163] + x[82] + x[63] + x[17] + 11 x[45] + x[61] + x[145]
   <= 11
 nod46: x[128] + x[94] + x[32] + x[65] + x[73] + x[100] + x[108] + x[123]
   + x[116] + x[105] + x[17] + 11 x[46] + x[151] + x[164] <= 11
 nod59: x[94] + x[65] + x[24] + x[141] + x[44] + x[146] + x[104] + x[117]
   + x[74] + x[95] + x[165] + x[38] + x[26] + x[110] + x[17] + 13 x[59]
   + x[168] + x[167] + x[157] <= 13
 nod132: x[128] + x[115] + x[131] + x[155] + x[134] + x[58] + x[100]
   + x[141] + x[119] + x[123] + x[66] + x[33] + x[124] + x[18] + 13 x[132]
   + x[173] + x[150] + x[47] <= 13
 nod174: x[23] + x[121] + x[73] + x[139] + x[76] + x[153] + x[141] + x[44]
   + x[22] + x[104] + x[42] + x[96] + x[91] + x[95] + x[30] + x[110]
   + x[18] + 12 x[174] + x[101] + x[145] + x[51] + x[162] <= 12
 nod138: x[128] + x[79] + x[169] + x[81] + x[119] + x[83] + x[113] + x[68]
   + x[120] + x[36] + x[26] + x[18] + x[82] + 12 x[138] + x[145] + x[31]
   + x[168] <= 12
 nod101: x[143] + x[65] + x[155] + x[73] + x[169] + x[92] + x[103] + x[160]
   + x[44] + x[119] + x[80] + x[96] + x[36] + x[26] + x[116] + x[21]
   + x[170] + x[19] + x[174] + 12 x[101] + x[137] + x[157] <= 12
 nod151: x[23] + x[155] + x[136] + x[57] + x[78] + x[97] + x[41] + x[71]
   + x[80] + x[68] + x[120] + x[20] + x[21] + x[19] + x[46] + 13 x[151]
   + x[99] + x[156] + x[168] <= 13
 nod173: x[77] + x[32] + x[44] + x[147] + x[123] + x[112] + x[20] + x[149]
   + x[154] + x[21] + x[132] + 10 x[173] + x[142] + x[55] <= 10
 nod99: x[73] + x[130] + x[24] + x[75] + x[50] + x[35] + x[36] + x[20]
   + x[151] + 10 x[99] + x[142] + x[43] + x[69] <= 10
 nod85: x[115] + x[87] + x[153] + x[41] + x[22] + x[29] + x[42] + x[53]
   + x[38] + x[120] + x[20] + x[163] + 9 x[85] + x[142] + x[61] <= 9
 nod86: x[67] + x[23] + x[135] + x[28] + x[37] + x[92] + x[44] + x[29]
   + x[50] + x[165] + x[68] + x[20] + x[21] + 10 x[86] + x[164] + x[168]
   <= 10
 nod150: x[67] + x[39] + x[89] + x[23] + x[94] + x[73] + x[62] + x[88]
   + x[81] + x[71] + x[107] + x[147] + x[112] + x[60] + x[21] + x[102]
   + x[132] + 12 x[150] + x[144] <= 12
 nod49: x[39] + x[23] + x[94] + x[169] + x[111] + x[62] + x[130] + x[171]
   + x[78] + x[41] + x[22] + x[96] + x[53] + x[95] + x[68] + x[54] + x[122]
   + 13 x[49] + x[162] <= 13
 nod142: x[73] + x[57] + x[22] + x[95] + x[133] + x[26] + x[110] + x[122]
   + x[173] + x[99] + x[85] + 8 x[142] <= 8
 nod61: x[65] + x[135] + x[72] + x[22] + x[29] + x[107] + x[53] + x[35]
   + x[26] + x[102] + x[45] + x[85] + 9 x[61] + x[51] + x[168] + x[47]
   <= 9
 nod144: x[129] + x[39] + x[23] + x[111] + x[109] + x[125] + x[126] + x[71]
   + x[42] + x[83] + x[147] + x[53] + x[95] + x[165] + x[52] + x[149]
   + x[150] + 10 x[144] <= 10
 nod145: x[67] + x[77] + x[23] + x[121] + x[134] + x[72] + x[76] + x[103]
   + x[109] + x[125] + x[80] + x[107] + x[53] + x[161] + x[66] + x[110]
   + x[45] + x[174] + x[138] + 14 x[145] + x[69] <= 14
 nod43: x[70] + x[23] + x[94] + x[65] + x[148] + x[57] + x[160] + x[166]
   + x[117] + x[113] + x[106] + x[120] + x[149] + x[154] + x[163] + x[90]
   + x[99] + 12 x[43] <= 12
 nod55: x[129] + x[67] + x[23] + x[34] + x[88] + x[100] + x[125] + x[160]
   + x[25] + x[83] + x[113] + x[165] + x[120] + x[60] + x[26] + x[173]
   + 12 x[55] + x[69] + x[159] + x[157] <= 12
 nod51: x[128] + x[79] + x[131] + x[72] + x[172] + x[24] + x[34] + x[97]
   + x[160] + x[41] + x[108] + x[44] + x[146] + x[152] + x[161] + x[36]
   + x[149] + x[170] + x[82] + x[174] + x[61] + 11 x[51] <= 11
 nod164: x[24] + x[160] + x[44] + x[25] + x[107] + x[75] + x[147] + x[68]
   + x[46] + x[86] + 8 x[164] + x[156] <= 8
 nod156: x[77] + x[65] + x[135] + x[155] + x[25] + x[29] + x[104] + x[75]
   + x[38] + x[120] + x[66] + x[151] + x[164] + 8 x[156] + x[47] <= 8
 nod137: x[118] + x[131] + x[58] + x[92] + x[100] + x[126] + x[81] + x[117]
   + x[133] + x[120] + x[26] + x[64] + x[101] + 9 x[137] <= 9
 nod31: x[67] + x[89] + x[135] + x[136] + x[58] + x[78] + x[48] + x[42]
   + x[91] + x[50] + x[114] + x[163] + x[138] + 10 x[31] + x[69] <= 10
 nod69: x[34] + x[175] + x[103] + x[146] + x[83] + x[91] + x[40] + x[99]
   + x[145] + x[55] + x[31] + 10 x[69] + x[168] + x[157] <= 10
 nod162: x[77] + x[32] + x[72] + x[73] + x[140] + x[103] + x[127] + x[117]
   + x[110] + x[122] + x[82] + x[174] + x[49] + 8 x[162] + x[159] <= 8
 nod159: x[32] + x[155] + x[130] + x[160] + x[80] + x[161] + x[68] + x[116]
   + x[55] + x[162] + 8 x[159] + x[168] <= 8
 nod168: x[70] + x[89] + x[172] + x[107] + x[84] + x[114] + x[33] + x[63]
   + x[59] + x[138] + x[151] + x[86] + x[61] + x[69] + x[159] + 12 x[168]
   <= 12
 nod167: x[39] + x[76] + x[37] + x[58] + x[92] + x[166] + x[50] + x[149]
   + x[59] + 7 x[167] <= 7
 nod47: x[70] + x[39] + x[135] + x[93] + x[158] + x[113] + x[161] + x[82]
   + x[132] + x[61] + x[156] + 7 x[47] + x[157] <= 7
 nod157: x[143] + x[94] + x[111] + x[62] + x[125] + x[41] + x[108] + x[44]
   + x[158] + x[96] + x[54] + x[60] + x[66] + x[116] + x[64] + x[82]
   + x[59] + x[101] + x[55] + x[69] + x[47] + 12 x[157] <= 12
Bounds
Binaries
 x[1] x[128] x[129] x[67] x[70] x[7] x[8] x[39] x[89] x[77] x[143] x[27]
 x[115] x[118] x[23] x[121] x[87] x[79] x[94] x[2] x[32] x[65] x[135]
 x[131] x[155] x[134] x[72] x[73] x[139] x[76] x[13] x[136] x[148] x[169]
 x[111] x[28] x[93] x[62] x[3] x[11] x[140] x[172] x[153] x[4] x[130] x[16]
 x[24] x[34] x[37] x[171] x[175] x[57] x[58] x[78] x[88] x[92] x[97] x[100]
 x[103] x[109] x[125] x[126] x[127] x[5] x[160] x[48] x[166] x[41] x[108]
 x[141] x[44] x[81] x[146] x[22] x[119] x[25] x[29] x[158] x[6] x[98] x[71]
 x[104] x[42] x[12] x[80] x[83] x[117] x[74] x[96] x[107] x[75] x[113]
 x[147] x[53] x[56] x[152] x[91] x[95] x[133] x[161] x[50] x[9] x[123]
 x[165] x[35] x[68] x[38] x[106] x[112] x[84] x[120] x[54] x[60] x[10]
 x[114] x[36] x[26] x[66] x[116] x[52] x[20] x[33] x[149] x[154] x[30]
 x[110] x[124] x[14] x[64] x[163] x[18] x[21] x[122] x[15] x[90] x[170]
 x[19] x[40] x[82] x[105] x[63] x[17] x[102] x[45] x[46] x[59] x[132]
 x[174] x[138] x[101] x[151] x[173] x[99] x[85] x[86] x[150] x[49] x[142]
 x[61] x[144] x[145] x[43] x[55] x[51] x[164] x[156] x[137] x[31] x[69]
 x[162] x[159] x[168] x[167] x[47] x[157]
End
