\ LP format - for model browsing. Use MPS format to capture full model detail.
Maximize
obj:   x[1] + x[130] + x[3] + x[5] + x[6] + x[136] + x[137] + x[14] + x[15]
   + x[147] + x[21] + x[27] + x[31] + x[33] + x[46] + x[49] + x[50] + x[57]
   + x[67] + x[68] + x[74] + x[76] + x[82] + x[91] + x[93] + x[97] + x[98]
   + x[99] + x[103] + x[107] + x[108] + x[109] + x[111] + x[112] + x[116]
   + x[124] + x[2] + x[131] + x[134] + x[138] + x[139] + x[145] + x[19]
   + x[24] + x[25] + x[35] + x[48] + x[54] + x[58] + x[61] + x[69] + x[70]
   + x[77] + x[79] + x[84] + x[86] + x[87] + x[117] + x[121] + x[122]
   + x[123] + x[125] + x[141] + x[148] + x[28] + x[29] + x[36] + x[37]
   + x[43] + x[53] + x[60] + x[62] + x[64] + x[65] + x[80] + x[83] + x[95]
   + x[105] + x[115] + x[4] + x[128] + x[129] + x[10] + x[143] + x[16]
   + x[18] + x[26] + x[55] + x[71] + x[81] + x[85] + x[92] + x[144]
   + x[104] + x[113] + x[114] + x[118] + x[146] + x[20] + x[149] + x[150]
   + x[39] + x[66] + x[72] + x[88] + x[106] + x[119] + x[133] + x[135]
   + x[8] + x[140] + x[22] + x[23] + x[32] + x[41] + x[44] + x[56] + x[75]
   + x[90] + x[127] + x[126] + x[7] + x[30] + x[94] + x[12] + x[40] + x[59]
   + x[63] + x[9] + x[142] + x[17] + x[42] + x[47] + x[52] + x[73] + x[100]
   + x[132] + x[11] + x[45] + x[102] + x[110] + x[89] + x[13] + x[34]
   + x[51] + x[78] + x[96] + x[101] + x[38] + x[120]
Subject To
 clq1: x[84] + x[53] + x[114] + x[56] <= 1
 clq2: x[46] + x[115] + x[114] + x[106] + x[22] <= 1
 clq3: x[77] + x[88] + x[140] + x[22] <= 1
 clq4: x[6] + x[21] + x[84] + x[41] <= 1
 clq5: x[82] + x[122] + x[83] + x[30] <= 1
 clq6: x[136] + x[80] + x[106] + x[142] <= 1
 clq7: x[53] + x[64] + x[4] + x[129] <= 1
 clq8: x[57] + x[66] + x[119] + x[59] <= 1
 clq9: x[103] + x[134] + x[4] + x[10] <= 1
 clq10: x[61] + x[83] + x[88] <= 1
 clq11: x[1] + x[14] + x[21] <= 1
 clq12: x[124] + x[84] + x[118] <= 1
 clq13: x[68] + x[125] + x[20] + x[30] <= 1
 clq14: x[1] + x[130] + x[136] + x[33] <= 1
 clq15: x[49] + x[69] + x[64] + x[42] <= 1
 clq16: x[111] + x[43] + x[88] + x[44] <= 1
 clq17: x[77] + x[29] + x[26] + x[75] <= 1
 clq18: x[84] + x[122] + x[119] + x[90] <= 1
 clq19: x[21] + x[53] + x[126] + x[94] <= 1
 clq20: x[31] + x[103] + x[61] + x[118] <= 1
 clq21: x[2] + x[134] + x[77] <= 1
 clq22: x[70] + x[141] + x[20] + x[102] <= 1
 clq23: x[5] + x[20] + x[39] + x[106] <= 1
 clq24: x[98] + x[16] + x[26] + x[22] <= 1
 clq25: x[4] + x[92] + x[104] + x[114] <= 1
 clq26: x[14] + x[53] + x[39] + x[127] + x[59] <= 1
 clq27: x[69] + x[121] + x[8] + x[7] + x[30] <= 1
 clq28: x[37] + x[104] + x[23] + x[100] <= 1
 clq29: x[130] + x[91] + x[54] <= 1
 clq30: x[136] + x[115] + x[150] <= 1
 clq31: x[137] + x[27] + x[134] <= 1
 clq32: x[49] + x[2] + x[25] + x[122] + x[125] <= 1
 clq33: x[111] + x[124] + x[58] + x[41] <= 1
 clq34: x[6] + x[19] + x[113] + x[75] <= 1
 clq35: x[86] + x[105] + x[81] + x[118] <= 1
 clq36: x[64] + x[113] + x[22] <= 1
 clq37: x[62] + x[20] + x[88] + x[142] <= 1
 clq38: x[31] + x[119] + x[52] <= 1
 clq39: x[1] + x[5] + x[15] + x[98] <= 1
 clq40: x[21] + x[2] + x[138] + x[48] <= 1
 clq41: x[103] + x[122] + x[44] <= 1
 clq42: x[109] + x[61] + x[42] + x[38] <= 1
 clq43: x[79] + x[115] + x[56] <= 1
 clq44: x[84] + x[113] + x[66] + x[94] <= 1
 clq45: x[136] + x[108] + x[141] <= 1
 clq46: x[74] + x[83] + x[51] <= 1
 clq47: x[134] + x[95] + x[72] <= 1
 clq48: x[68] + x[106] + x[9] <= 1
 clq49: x[130] + x[14] + x[73] + x[45] <= 1
 clq50: x[27] + x[54] + x[10] + x[100] <= 1
 clq51: x[50] + x[16] + x[17] <= 1
 clq52: x[137] + x[97] + x[79] + x[37] <= 1
 clq53: x[99] + x[116] + x[77] <= 1
 clq54: x[64] + x[114] + x[44] <= 1
 clq55: x[54] + x[36] + x[55] + x[102] <= 1
 clq56: x[46] + x[117] + x[26] + x[150] <= 1
 clq57: x[50] + x[119] + x[75] + x[94] <= 1
 clq58: x[31] + x[106] + x[12] + x[73] <= 1
 clq59: x[103] + x[48] + x[58] + x[56] <= 1
 clq60: x[111] + x[24] + x[141] <= 1
 clq61: x[68] + x[124] + x[4] <= 1
 clq62: x[131] + x[43] + x[59] <= 1
 clq63: x[97] + x[19] + x[71] + x[52] <= 1
 clq64: x[70] + x[62] + x[65] <= 1
 clq65: x[86] + x[18] + x[127] <= 1
 clq66: x[29] + x[85] + x[41] + x[7] <= 1
 clq67: x[109] + x[77] + x[83] + x[135] <= 1
 clq68: x[33] + x[124] + x[95] + x[22] <= 1
 clq69: x[15] + x[70] + x[105] <= 1
 clq70: x[139] + x[81] + x[142] + x[45] <= 1
 clq71: x[114] + x[72] + x[88] + x[96] <= 1
 clq72: x[136] + x[98] + x[23] <= 1
 clq73: x[108] + x[30] + x[13] <= 1
 clq74: x[1] + x[147] + x[57] <= 1
 clq75: x[21] + x[145] + x[144] <= 1
 clq76: x[49] + x[40] + x[89] <= 1
 clq77: x[99] + x[16] + x[9] <= 1
 clq78: x[116] + x[24] + x[84] + x[71] <= 1
 clq79: x[107] + x[134] + x[115] <= 1
 clq80: x[15] + x[131] + x[25] + x[53] <= 1
 clq81: x[3] + x[61] + x[60] + x[65] <= 1
 clq82: x[130] + x[74] + x[64] <= 1
 clq83: x[27] + x[141] + x[51] <= 1
 clq84: x[5] + x[49] + x[95] + x[18] <= 1
 clq85: x[67] + x[22] + x[23] + x[63] <= 1
 clq86: x[136] + x[91] + x[43] + x[89] <= 1
 clq87: x[93] + x[56] + x[30] <= 1
 clq88: x[111] + x[22] + x[9] + x[100] <= 1
 clq89: x[86] + x[59] + x[110] <= 1
 clq90: x[77] + x[105] + x[4] <= 1
 clq91: x[83] + x[4] + x[128] + x[55] <= 1
 clq92: x[87] + x[118] + x[44] + x[42] <= 1
 clq93: x[148] + x[150] + x[90] <= 1
 clq94: x[108] + x[119] + x[102] <= 1
 clq95: x[68] + x[112] + x[127] + x[132] <= 1
 clq96: x[109] + x[134] + x[126] <= 1
 clq97: x[124] + x[122] + x[16] + x[13] <= 1
 clq98: x[140] + x[90] + x[17] + x[38] <= 1
 clq99: x[1] + x[6] + x[137] <= 1
 clq100: x[147] + x[21] + x[106] <= 1
 clq101: x[31] + x[10] + x[149] <= 1
 clq102: x[76] + x[99] + x[36] + x[12] <= 1
 clq103: x[103] + x[115] + x[149] + x[51] <= 1
 clq104: x[21] + x[19] + x[66] + x[110] <= 1
 clq105: x[50] + x[69] + x[148] + x[26] + x[81] <= 1
 clq106: x[79] + x[41] + x[126] <= 1
 clq107: x[3] + x[97] + x[29] + x[62] <= 1
 clq108: x[3] + x[98] + x[37] + x[64] <= 1
 clq109: x[35] + x[54] + x[88] + x[75] <= 1
 clq110: x[64] + x[133] + x[52] <= 1
 clq111: x[74] + x[75] + x[142] + x[96] <= 1
 clq112: x[130] + x[82] + x[135] + x[101] <= 1
 clq113: x[14] + x[116] + x[138] <= 1
 clq114: x[27] + x[145] + x[86] <= 1
 clq115: x[57] + x[48] + x[114] + x[7] <= 1
 clq116: x[46] + x[93] + x[144] <= 1
 clq117: x[138] + x[25] + x[84] + x[65] <= 1
 clq118: x[61] + x[141] + x[113] <= 1
 clq119: x[33] + x[70] + x[77] + x[150] <= 1
 clq120: x[147] + x[122] + x[105] <= 1
 clq121: x[104] + x[30] + x[110] <= 1
 clq122: x[131] + x[39] + x[94] <= 1
 clq123: x[5] + x[104] + x[72] + x[119] <= 1
 clq124: x[44] + x[40] + x[59] + x[9] <= 1
 clq125: x[10] + x[106] + x[102] <= 1
 clq126: x[1] + x[31] + x[68] <= 1
 clq127: x[137] + x[146] + x[32] <= 1
 clq128: x[67] + x[116] + x[18] <= 1
 clq129: x[91] + x[107] + x[133] + x[41] <= 1
 clq130: x[108] + x[43] + x[56] <= 1
 clq131: x[111] + x[134] + x[37] + x[92] <= 1
 clq132: x[2] + x[61] + x[79] <= 1
 clq133: x[117] + x[88] + x[45] <= 1
 clq134: x[125] + x[100] + x[38] <= 1
 clq135: x[77] + x[53] + x[66] + x[52] <= 1
 clq136: x[60] + x[128] + x[51] <= 1
 clq137: x[15] + x[55] + x[13] <= 1
 clq138: x[6] + x[118] + x[140] <= 1
 clq139: x[115] + x[23] + x[12] <= 1
 clq140: x[136] + x[99] + x[127] <= 1
 clq141: x[27] + x[103] + x[12] + x[132] <= 1
 clq142: x[14] + x[109] + x[34] <= 1
 clq143: x[49] + x[86] + x[135] <= 1
 clq144: x[93] + x[71] + x[17] <= 1
 clq145: x[2] + x[145] + x[24] <= 1
 clq146: x[33] + x[19] + x[141] <= 1
 clq147: x[136] + x[123] + x[36] + x[42] <= 1
 clq148: x[29] + x[143] + x[102] <= 1
 clq149: x[62] + x[10] + x[126] <= 1
 clq150: x[83] + x[144] + x[100] <= 1
 clq151: x[43] + x[95] + x[129] <= 1
 clq152: x[54] + x[105] + x[32] <= 1
 clq153: x[25] + x[16] + x[78] <= 1
 clq154: x[85] + x[30] + x[45] <= 1
 clq155: x[39] + x[22] + x[30] + x[142] <= 1
 clq156: x[94] + x[89] + x[13] <= 1
 clq157: x[113] + x[20] + x[96] <= 1
 clq158: x[1] + x[46] + x[99] <= 1
 clq159: x[130] + x[26] + x[146] <= 1
 clq160: x[137] + x[147] + x[84] <= 1
 clq161: x[31] + x[60] + x[7] <= 1
 clq162: x[57] + x[4] + x[18] <= 1
 clq163: x[67] + x[59] + x[78] <= 1
 clq164: x[74] + x[43] + x[55] <= 1
 clq165: x[76] + x[146] + x[23] <= 1
 clq166: x[82] + x[131] + x[69] <= 1
 clq167: x[108] + x[72] + x[73] <= 1
 clq168: x[48] + x[141] <= 1
 clq169: x[2] + x[70] + x[117] <= 1
 clq170: x[79] + x[40] + x[101] <= 1
 clq171: x[6] + x[121] + x[56] <= 1
 clq172: x[14] + x[125] + x[106] <= 1
 clq173: x[128] + x[119] + x[8] <= 1
 clq174: x[122] + x[81] + x[132] <= 1
 clq175: x[114] + x[78] + x[38] <= 1
 clq176: x[109] + x[143] + x[118] <= 1
 clq177: x[134] + x[150] + x[34] <= 1
 clq178: x[15] + x[60] + x[26] + x[90] <= 1
 clq179: x[115] + x[127] + x[17] <= 1
 clq180: x[5] + x[19] + x[25] + x[146] <= 1
 clq181: x[21] + x[29] + x[95] <= 1
 clq182: x[49] + x[50] + x[35] <= 1
 clq183: x[68] + x[18] + x[144] + x[102] <= 1
 clq184: x[76] + x[91] + x[140] <= 1
 clq185: x[93] + x[48] + x[53] <= 1
 clq186: x[97] + x[24] + x[88] <= 1
 clq187: x[49] + x[98] + x[105] <= 1
 clq188: x[107] + x[13] + x[51] <= 1
 clq189: x[111] + x[83] + x[42] <= 1
 clq190: x[108] + x[116] + x[16] <= 1
 clq191: x[124] + x[65] + x[96] <= 1
 clq192: x[145] + x[8] + x[41] <= 1
 clq193: x[61] + x[92] + x[94] <= 1
 clq194: x[21] + x[86] + x[52] <= 1
 clq195: x[62] + x[41] + x[44] <= 1
 clq196: x[123] + x[64] + x[7] <= 1
 clq197: x[80] + x[89] + x[34] <= 1
 clq198: x[117] + x[71] + x[63] <= 1
 clq199: x[14] + x[4] + x[81] + x[85] <= 1
 clq200: x[6] + x[133] + x[126] <= 1
 clq201: x[1] + x[27] + x[109] <= 1
 clq202: x[130] + x[113] + x[9] <= 1
 clq203: x[137] + x[138] + x[133] <= 1
 clq204: x[147] + x[131] + x[88] <= 1
 clq205: x[31] + x[129] <= 1
 clq206: x[46] + x[98] + x[55] + x[75] <= 1
 clq207: x[57] + x[121] + x[20] <= 1
 clq208: x[67] + x[122] + x[135] <= 1
 clq209: x[74] + x[141] + x[127] <= 1
 clq210: x[82] + x[103] + x[45] <= 1
 clq211: x[139] + x[64] + x[32] <= 1
 clq212: x[27] + x[69] + x[104] + x[40] <= 1
 clq213: x[70] + x[125] + x[39] <= 1
 clq214: x[46] + x[79] + x[110] <= 1
 clq215: x[145] + x[84] + x[72] <= 1
 clq216: x[50] + x[37] + x[38] <= 1
 clq217: x[62] + x[128] + x[47] <= 1
 clq218: x[5] + x[10] + x[150] <= 1
 clq219: x[116] + x[139] + x[114] + x[119] <= 1
 clq220: x[118] + x[110] + x[101] <= 1
 clq221: x[131] + x[66] + x[132] <= 1
 clq222: x[105] + x[129] + x[106] <= 1
 clq223: x[48] + x[55] + x[23] <= 1
 clq224: x[56] + x[100] + x[45] <= 1
 clq225: x[29] + x[126] + x[59] <= 1
 clq226: x[31] + x[128] + x[71] <= 1
 clq227: x[33] + x[99] + x[2] <= 1
 clq228: x[91] + x[65] <= 1
 clq229: x[93] + x[140] + x[7] <= 1
 clq230: x[97] + x[83] + x[94] <= 1
 clq231: x[107] + x[58] + x[143] <= 1
 clq232: x[111] + x[36] <= 1
 clq233: x[109] + x[112] + x[40] + x[51] <= 1
 clq234: x[124] + x[53] <= 1
 clq235: x[134] + x[44] + x[101] <= 1
 clq236: x[24] + x[35] + x[45] <= 1
 clq237: x[54] + x[148] <= 1
 clq238: x[61] + x[53] + x[115] <= 1
 clq239: x[15] + x[86] + x[17] <= 1
 clq240: x[80] + x[16] + x[72] <= 1
 clq241: x[136] + x[95] + x[52] <= 1
 clq242: x[67] + x[44] + x[90] + x[52] <= 1
 clq243: x[68] + x[91] + x[63] <= 1
 clq244: x[9] + x[142] + x[73] <= 1
 clq245: x[130] + x[124] + x[149] + x[89] <= 1
 clq246: x[108] + x[96] <= 1
 clq247: x[1] + x[3] + x[49] <= 1
 clq248: x[137] + x[85] + x[150] <= 1
 clq249: x[147] + x[79] + x[102] <= 1
 clq250: x[57] + x[60] + x[38] <= 1
 clq251: x[74] + x[16] + x[66] <= 1
 clq252: x[76] + x[86] + x[34] <= 1
 clq253: x[82] + x[107] + x[84] <= 1
 clq254: x[103] + x[65] + x[80] <= 1
 clq255: x[99] + x[111] + x[39] <= 1
 clq256: x[138] + x[106] + x[41] <= 1
 clq257: x[19] + x[36] + x[72] <= 1
 clq258: x[54] + x[43] + x[114] + x[132] <= 1
 clq259: x[111] + x[70] + x[144] <= 1
 clq260: x[87] + x[23] + x[120] <= 1
 clq261: x[117] + x[20] + x[47] <= 1
 clq262: x[121] + x[37] <= 1
 clq263: x[2] + x[35] + x[123] <= 1
 clq264: x[77] + x[125] + x[59] <= 1
 clq265: x[141] + x[115] + x[42] <= 1
 clq266: x[48] + x[148] + x[22] + x[12] <= 1
 clq267: x[28] + x[10] + x[92] + x[133] <= 1
 clq268: x[33] + x[93] + x[118] <= 1
 clq269: x[29] + x[135] + x[100] <= 1
 clq270: x[32] + x[127] + x[13] <= 1
 clq271: x[104] + x[20] + x[56] <= 1
 clq272: x[30] + x[73] <= 1
 clq273: x[134] + x[94] + x[96] <= 1
 clq274: x[6] + x[4] <= 1
 clq275: x[21] + x[37] <= 1
 clq276: x[31] + x[117] + x[75] <= 1
 clq277: x[1] + x[50] + x[97] <= 1
 clq278: x[21] + x[108] + x[113] <= 1
 clq279: x[145] + x[129] <= 1
 clq280: x[24] + x[79] + x[119] <= 1
 clq281: x[25] + x[95] + x[128] <= 1
 clq282: x[58] + x[61] + x[87] <= 1
 clq283: x[69] + x[83] + x[142] <= 1
 clq284: x[121] + x[105] <= 1
 clq285: x[91] + x[122] + x[71] <= 1
 clq286: x[82] + x[62] + x[129] <= 1
 clq287: x[64] + x[143] + x[39] <= 1
 clq288: x[50] + x[18] + x[126] <= 1
 clq289: x[4] + x[26] + x[118] <= 1
 clq290: x[37] + x[81] + x[63] <= 1
 clq291: x[5] + x[66] + x[88] <= 1
 clq292: x[35] + x[65] + x[8] + x[12] <= 1
 clq293: x[70] + x[80] + x[140] <= 1
 clq294: x[57] + x[30] + x[78] <= 1
 clq295: x[98] + x[131] + x[17] <= 1
 clq296: x[3] + x[36] + x[53] <= 1
 clq297: x[6] + x[135] + x[44] <= 1
 clq298: x[136] + x[7] <= 1
 clq299: x[137] + x[67] + x[101] <= 1
 clq300: x[14] + x[100] <= 1
 clq301: x[15] + x[42] + x[110] <= 1
 clq302: x[147] + x[103] + x[73] <= 1
 clq303: x[68] + x[105] + x[143] <= 1
 clq304: x[74] + x[100] + x[102] <= 1
 clq305: x[76] + x[120] <= 1
 clq306: x[112] + x[95] + x[38] <= 1
 clq307: x[116] + x[59] + x[13] <= 1
 clq308: x[138] + x[34] <= 1
 clq309: x[131] + x[145] + x[89] <= 1
 clq310: x[19] + x[144] <= 1
 clq311: x[58] + x[84] + x[123] <= 1
 clq312: x[121] + x[10] + x[11] <= 1
 clq313: x[125] + x[64] + x[96] <= 1
 clq314: x[97] + x[141] + x[41] <= 1
 clq315: x[85] + x[127] + x[94] <= 1
 clq316: x[104] + x[150] + x[9] <= 1
 clq317: x[33] + x[92] + x[106] <= 1
 clq318: x[6] + x[23] + x[32] <= 1
 clq319: x[10] + x[56] + x[47] <= 1
 clq320: x[90] + x[45] + x[51] <= 1
 clq321: x[130] + x[34] <= 1
 clq322: x[136] + x[139] <= 1
 clq323: x[14] + x[122] + x[26] <= 1
 clq324: x[27] + x[108] + x[133] <= 1
 clq325: x[31] + x[144] <= 1
 clq326: x[46] + x[77] <= 1
 clq327: x[46] + x[49] + x[142] <= 1
 clq328: x[76] + x[79] + x[81] <= 1
 clq329: x[3] + x[93] + x[83] <= 1
 clq330: x[99] + x[18] <= 1
 clq331: x[107] + x[43] + x[26] <= 1
 clq332: x[1] + x[124] <= 1
 clq333: x[82] + x[134] + x[7] <= 1
 clq334: x[111] + x[138] + x[140] <= 1
 clq335: x[19] + x[16] + x[92] <= 1
 clq336: x[24] + x[86] + x[148] <= 1
 clq337: x[25] + x[54] <= 1
 clq338: x[61] + x[119] + x[40] <= 1
 clq339: x[69] + x[29] + x[16] <= 1
 clq340: x[60] + x[146] + x[44] <= 1
 clq341: x[62] + x[149] <= 1
 clq342: x[31] + x[115] + x[135] <= 1
 clq343: x[129] + x[55] + x[81] <= 1
 clq344: x[71] + x[102] + x[120] <= 1
 clq345: x[50] + x[113] + x[78] <= 1
 clq346: x[49] + x[88] + x[63] <= 1
 clq347: x[75] + x[17] <= 1
 clq348: x[99] + x[126] <= 1
 clq349: x[98] + x[30] <= 1
 clq350: x[130] + x[66] <= 1
 clq351: x[5] + x[37] <= 1
 clq352: x[136] + x[11] <= 1
 clq353: x[137] + x[145] + x[110] <= 1
 clq354: x[15] + x[2] <= 1
 clq355: x[147] + x[72] <= 1
 clq356: x[21] + x[36] + x[78] <= 1
 clq357: x[57] + x[116] + x[70] <= 1
 clq358: x[68] + x[150] + x[101] <= 1
 clq359: x[74] + x[73] <= 1
 clq360: x[91] + x[125] + x[56] <= 1
 clq361: x[103] + x[20] + x[32] <= 1
 clq362: x[109] + x[106] <= 1
 clq363: x[112] + x[85] <= 1
 clq364: x[124] + x[8] + x[59] <= 1
 clq365: x[139] + x[9] <= 1
 clq366: x[82] + x[25] + x[52] <= 1
 clq367: x[48] + x[39] <= 1
 clq368: x[97] + x[54] + x[34] <= 1
 clq369: x[77] + x[141] <= 1
 clq370: x[46] + x[84] + x[89] <= 1
 clq371: x[87] + x[105] <= 1
 clq372: x[117] + x[114] <= 1
 clq373: x[121] + x[132] <= 1
 clq374: x[123] + x[9] + x[47] <= 1
 clq375: x[53] + x[60] + x[41] <= 1
 clq376: x[62] + x[38] <= 1
 clq377: x[138] + x[4] + x[144] <= 1
 clq378: x[33] + x[128] <= 1
 clq379: x[134] + x[18] <= 1
 clq380: x[35] + x[118] + x[17] <= 1
 clq381: x[71] + x[149] + x[22] <= 1
 clq382: x[72] + x[23] <= 1
 clq383: x[75] + x[51] <= 1
 clq384: x[123] + x[105] + x[90] <= 1
 clq385: x[127] + x[96] <= 1
 clq386: x[69] + x[87] + x[126] <= 1
 clq387: x[30] + x[42] <= 1
 clq388: x[136] + x[57] + x[13] <= 1
 clq389: x[1] + x[67] + x[76] <= 1
 clq390: x[130] + x[11] <= 1
 clq391: x[3] + x[5] + x[43] <= 1
 clq392: x[14] + x[31] + x[98] <= 1
 clq393: x[15] + x[39] + x[126] + x[100] <= 1
 clq394: x[147] + x[128] + x[94] <= 1
 clq395: x[27] + x[42] <= 1
 clq396: x[74] + x[62] + x[7] <= 1
 clq397: x[93] + x[2] + x[131] <= 1
 clq398: x[99] + x[122] <= 1
 clq399: x[107] + x[150] + x[66] <= 1
 clq400: x[108] + x[77] <= 1
 clq401: x[109] + x[80] <= 1
 clq402: x[112] + x[115] <= 1
 clq403: x[139] + x[69] + x[106] <= 1
 clq404: x[19] + x[89] + x[38] <= 1
 clq405: x[19] + x[24] + x[104] <= 1
 clq406: x[48] + x[51] + x[96] <= 1
 clq407: x[58] + x[65] + x[47] <= 1
 clq408: x[61] + x[143] <= 1
 clq409: x[33] + x[86] <= 1
 clq410: x[117] + x[55] <= 1
 clq411: x[134] + x[121] + x[140] <= 1
 clq412: x[141] + x[8] + x[63] <= 1
 clq413: x[148] + x[73] <= 1
 clq414: x[107] + x[29] + x[113] <= 1
 clq415: x[37] + x[146] <= 1
 clq416: x[64] + x[132] <= 1
 clq417: x[31] + x[83] + x[34] <= 1
 clq418: x[99] + x[95] + x[56] <= 1
 clq419: x[18] + x[114] <= 1
 clq420: x[27] + x[85] + x[88] <= 1
 clq421: x[129] + x[119] <= 1
 clq422: x[133] + x[23] + x[101] <= 1
 clq423: x[36] + x[75] <= 1
 clq424: x[6] + x[70] + x[127] <= 1
 clq425: x[109] + x[30] <= 1
 clq426: x[98] + x[142] + x[120] <= 1
 clq427: x[130] + x[50] + x[143] <= 1
 clq428: x[137] + x[80] <= 1
 clq429: x[21] + x[42] <= 1
 clq430: x[68] + x[108] + x[100] <= 1
 clq431: x[1] + x[91] + x[116] <= 1
 clq432: x[103] + x[102] <= 1
 clq433: x[111] + x[12] <= 1
 clq434: x[112] + x[119] + x[45] <= 1
 clq435: x[124] + x[102] <= 1
 clq436: x[145] + x[37] + x[65] <= 1
 clq437: x[25] + x[73] <= 1
 clq438: x[145] + x[54] + x[18] <= 1
 clq439: x[61] + x[72] <= 1
 clq440: x[93] + x[77] + x[11] <= 1
 clq441: x[79] + x[148] + x[4] <= 1
 clq442: x[84] + x[55] + x[135] <= 1
 clq443: x[86] + x[122] <= 1
 clq444: x[117] + x[52] <= 1
 clq445: x[125] + x[144] + x[40] <= 1
 clq446: x[108] + x[28] + x[90] <= 1
 clq447: x[53] + x[16] <= 1
 clq448: x[64] + x[17] <= 1
 clq449: x[115] + x[32] + x[41] <= 1
 clq450: x[114] + x[75] <= 1
 clq451: x[118] + x[146] <= 1
 clq452: x[20] + x[132] <= 1
 clq453: x[149] + x[150] + x[30] <= 1
 clq454: x[86] + x[129] + x[22] <= 1
 clq455: x[59] + x[38] <= 1
 clq456: x[73] + x[110] <= 1
 clq457: x[106] + x[13] <= 1
 clq458: x[5] + x[97] + x[125] <= 1
 clq459: x[136] + x[111] + x[17] <= 1
 clq460: x[137] + x[49] <= 1
 clq461: x[14] + x[67] <= 1
 clq462: x[147] + x[138] <= 1
 clq463: x[21] + x[104] <= 1
 clq464: x[33] + x[103] <= 1
 clq465: x[74] + x[76] <= 1
 clq466: x[109] + x[66] <= 1
 clq467: x[124] + x[43] + x[92] <= 1
 clq468: x[2] + x[139] <= 1
 clq469: x[24] + x[80] <= 1
 clq470: x[25] + x[10] <= 1
 clq471: x[48] + x[13] <= 1
 clq472: x[58] + x[141] <= 1
 clq473: x[61] + x[105] <= 1
 clq474: x[87] + x[96] <= 1
 clq475: x[117] + x[81] + x[59] <= 1
 clq476: x[121] + x[12] <= 1
 clq477: x[122] + x[72] <= 1
 clq478: x[29] + x[18] + x[9] <= 1
 clq479: x[36] + x[9] + x[52] <= 1
 clq480: x[108] + x[53] + x[83] <= 1
 clq481: x[60] + x[85] + x[110] <= 1
 clq482: x[62] + x[118] <= 1
 clq483: x[64] + x[38] <= 1
 clq484: x[95] + x[42] + x[45] <= 1
 clq485: x[130] + x[15] + x[128] <= 1
 clq486: x[26] + x[114] <= 1
 clq487: x[71] + x[23] <= 1
 clq488: x[113] + x[63] <= 1
 clq489: x[146] + x[7] <= 1
 clq490: x[131] + x[20] <= 1
 clq491: x[88] + x[132] <= 1
 clq492: x[122] + x[140] + x[51] <= 1
 clq493: x[67] + x[75] <= 1
 clq494: x[127] + x[7] <= 1
 clq495: x[94] + x[132] <= 1
 clq496: x[35] + x[142] <= 1
 clq497: x[34] + x[78] <= 1
 clq498: x[139] + x[66] + x[120] <= 1
 clq499: x[3] + x[80] + x[105] <= 1
 clq500: x[6] + x[8] <= 1
 clq501: x[137] + x[113] <= 1
 clq502: x[14] + x[37] + x[62] <= 1
 clq503: x[147] + x[118] <= 1
 clq504: x[21] + x[142] <= 1
 clq505: x[27] + x[106] + x[78] <= 1
 clq506: x[33] + x[135] <= 1
 clq507: x[46] + x[103] <= 1
 clq508: x[49] + x[4] <= 1
 clq509: x[57] + x[71] <= 1
 clq510: x[68] + x[47] + x[13] <= 1
 clq511: x[74] + x[134] + x[87] <= 1
 clq512: x[76] + x[70] <= 1
 clq513: x[82] + x[75] <= 1
 clq514: x[91] + x[94] <= 1
 clq515: x[99] + x[26] + x[11] <= 1
 clq516: x[109] + x[138] + x[86] <= 1
 clq517: x[112] + x[20] <= 1
 clq518: x[116] + x[141] + x[26] <= 1
 clq519: x[2] + x[84] <= 1
 clq520: x[24] + x[22] <= 1
 clq521: x[50] + x[25] + x[55] <= 1
 clq522: x[48] + x[28] + x[102] <= 1
 clq523: x[54] + x[30] <= 1
 clq524: x[58] + x[40] <= 1
 clq525: x[137] + x[61] <= 1
 clq526: x[77] + x[115] <= 1
 clq527: x[79] + x[146] <= 1
 clq528: x[121] + x[104] <= 1
 clq529: x[123] + x[88] <= 1
 clq530: x[130] + x[148] + x[71] <= 1
 clq531: x[64] + x[73] <= 1
 clq532: x[97] + x[10] <= 1
 clq533: x[143] + x[90] <= 1
 clq534: x[16] + x[101] <= 1
 clq535: x[144] + x[114] <= 1
 clq536: x[149] + x[119] <= 1
 clq537: x[39] + x[32] <= 1
 clq538: x[15] + x[72] + x[56] <= 1
 clq539: x[8] + x[23] <= 1
 clq540: x[44] + x[34] <= 1
 clq541: x[127] + x[12] <= 1
 clq542: x[19] + x[63] <= 1
 clq543: x[89] + x[101] <= 1
 clq544: x[116] + x[132] + x[96] <= 1
 clq545: x[1] + x[103] <= 1
 clq546: x[5] + x[70] + x[83] <= 1
 clq547: x[6] + x[136] <= 1
 clq548: x[147] + x[41] <= 1
 clq549: x[21] + x[38] <= 1
 clq550: x[33] + x[82] <= 1
 clq551: x[46] + x[123] <= 1
 clq552: x[49] + x[110] <= 1
 clq553: x[57] + x[42] <= 1
 clq554: x[67] + x[92] <= 1
 clq555: x[76] + x[44] <= 1
 clq556: x[91] + x[95] + x[59] <= 1
 clq557: x[93] + x[87] + x[39] <= 1
 clq558: x[98] + x[4] <= 1
 clq559: x[107] + x[32] + x[30] <= 1
 clq560: x[111] + x[120] <= 1
 clq561: x[112] + x[128] <= 1
 clq562: x[124] + x[104] <= 1
 clq563: x[2] + x[54] + x[69] <= 1
 clq564: x[131] + x[35] + x[102] <= 1
 clq565: x[24] + x[53] + x[135] <= 1
 clq566: x[58] + x[146] <= 1
 clq567: x[61] + x[29] <= 1
 clq568: x[77] + x[7] <= 1
 clq569: x[79] + x[129] <= 1
 clq570: x[84] + x[140] <= 1
 clq571: x[117] + x[90] <= 1
 clq572: x[121] + x[122] + x[114] <= 1
 clq573: x[36] + x[118] <= 1
 clq574: x[43] + x[133] <= 1
 clq575: x[60] + x[81] + x[119] <= 1
 clq576: x[64] + x[92] + x[149] <= 1
 clq577: x[65] + x[40] <= 1
 clq578: x[115] + x[73] <= 1
 clq579: x[48] + x[10] + x[85] <= 1
 clq580: x[143] + x[75] <= 1
 clq581: x[16] + x[12] <= 1
 clq582: x[144] + x[66] <= 1
 clq583: x[113] + x[126] <= 1
 clq584: x[20] + x[17] <= 1
 clq585: x[88] + x[59] + x[51] <= 1
 clq586: x[65] + x[22] + x[11] <= 1
 clq587: x[138] + x[23] + x[94] <= 1
 clq588: x[67] + x[80] + x[127] <= 1
 clq589: x[99] + x[105] + x[63] <= 1
 clq590: x[129] + x[142] + x[17] <= 1
 clq591: x[68] + x[52] <= 1
 clq592: x[106] + x[89] <= 1
 clq593: x[1] + x[74] + x[93] <= 1
 clq594: x[3] + x[141] <= 1
 clq595: x[6] + x[46] <= 1
 clq596: x[136] + x[37] <= 1
 clq597: x[137] + x[124] + x[56] <= 1
 clq598: x[14] + x[134] <= 1
 clq599: x[147] + x[111] <= 1
 clq600: x[21] + x[8] <= 1
 clq601: x[27] + x[36] <= 1
 clq602: x[31] + x[13] <= 1
 clq603: x[33] + x[38] <= 1
 clq604: x[49] + x[111] + x[45] <= 1
 clq605: x[57] + x[12] <= 1
 clq606: x[76] + x[103] + x[25] <= 1
 clq607: x[82] + x[98] <= 1
 clq608: x[97] + x[109] <= 1
 clq609: x[112] + x[126] <= 1
 clq610: x[139] + x[77] <= 1
 clq611: x[145] + x[34] + x[38] <= 1
 clq612: x[2] + x[19] <= 1
 clq613: x[58] + x[4] <= 1
 clq614: x[61] + x[100] <= 1
 clq615: x[79] + x[113] <= 1
 clq616: x[84] + x[16] <= 1
 clq617: x[117] + x[101] <= 1
 clq618: x[123] + x[62] <= 1
 clq619: x[33] + x[125] <= 1
 clq620: x[148] + x[42] <= 1
 clq621: x[3] + x[28] + x[95] <= 1
 clq622: x[29] + x[110] <= 1
 clq623: x[43] + x[75] <= 1
 clq624: x[115] + x[72] <= 1
 clq625: x[128] + x[143] <= 1
 clq626: x[138] + x[18] + x[7] <= 1
 clq627: x[5] + x[55] <= 1
 clq628: x[144] + x[41] <= 1
 clq629: x[104] + x[44] <= 1
 clq630: x[118] + x[133] + x[51] <= 1
 clq631: x[68] + x[109] + x[146] <= 1
 clq632: x[77] + x[37] + x[20] <= 1
 clq633: x[150] + x[135] <= 1
 clq634: x[140] + x[40] <= 1
 clq635: x[36] + x[90] <= 1
 clq636: x[114] + x[47] + x[73] <= 1
 clq637: x[52] + x[100] <= 1
 clq638: x[71] + x[89] <= 1
 clq639: x[78] + x[96] <= 1
 clq640: x[132] + x[120] <= 1
 clq641: x[6] + x[82] <= 1
 clq642: x[136] + x[131] <= 1
 clq643: x[14] + x[78] <= 1
 clq644: x[147] + x[81] <= 1
 clq645: x[21] + x[60] <= 1
 clq646: x[27] + x[55] + x[8] <= 1
 clq647: x[31] + x[113] <= 1
 clq648: x[46] + x[44] <= 1
 clq649: x[50] + x[23] <= 1
 clq650: x[57] + x[84] <= 1
 clq651: x[91] + x[19] <= 1
 clq652: x[97] + x[98] + x[117] <= 1
 clq653: x[107] + x[62] <= 1
 clq654: x[112] + x[66] <= 1
 clq655: x[134] + x[75] <= 1
 clq656: x[124] + x[139] + x[70] <= 1
 clq657: x[99] + x[24] + x[120] <= 1
 clq658: x[35] + x[144] <= 1
 clq659: x[54] + x[85] <= 1
 clq660: x[58] + x[90] <= 1
 clq661: x[61] + x[150] <= 1
 clq662: x[69] + x[22] <= 1
 clq663: x[79] + x[52] + x[13] <= 1
 clq664: x[86] + x[101] <= 1
 clq665: x[121] + x[41] <= 1
 clq666: x[123] + x[106] <= 1
 clq667: x[125] + x[119] <= 1
 clq668: x[141] + x[140] <= 1
 clq669: x[148] + x[88] <= 1
 clq670: x[147] + x[29] + x[64] <= 1
 clq671: x[43] + x[115] <= 1
 clq672: x[53] + x[63] <= 1
 clq673: x[83] + x[104] <= 1
 clq674: x[4] + x[143] <= 1
 clq675: x[128] + x[39] <= 1
 clq676: x[10] + x[110] <= 1
 clq677: x[16] + x[149] + x[56] <= 1
 clq678: x[72] + x[42] <= 1
 clq679: x[32] + x[100] <= 1
 clq680: x[104] + x[127] <= 1
 clq681: x[126] + x[12] <= 1
 clq682: x[30] + x[89] <= 1
 clq683: x[7] + x[94] <= 1
 clq684: x[40] + x[132] <= 1
 clq685: x[41] + x[9] <= 1
 clq686: x[130] + x[133] + x[142] <= 1
 clq687: x[15] + x[96] <= 1
 clq688: x[1] + x[107] <= 1
 clq689: x[5] + x[108] <= 1
 clq690: x[6] + x[22] <= 1
 clq691: x[136] + x[83] <= 1
 clq692: x[137] + x[69] <= 1
 clq693: x[14] + x[47] <= 1
 clq694: x[21] + x[43] <= 1
 clq695: x[31] + x[2] <= 1
 clq696: x[33] + x[115] <= 1
 clq697: x[46] + x[148] <= 1
 clq698: x[49] + x[48] <= 1
 clq699: x[50] + x[67] <= 1
 clq700: x[57] + x[118] <= 1
 clq701: x[74] + x[94] + x[12] <= 1
 clq702: x[76] + x[61] <= 1
 clq703: x[82] + x[146] <= 1
 clq704: x[91] + x[26] <= 1
 clq705: x[103] + x[140] <= 1
 clq706: x[1] + x[108] + x[112] <= 1
 clq707: x[116] + x[17] <= 1
 clq708: x[131] + x[42] <= 1
 clq709: x[134] + x[135] <= 1
 clq710: x[145] + x[73] <= 1
 clq711: x[19] + x[129] <= 1
 clq712: x[25] + x[63] <= 1
 clq713: x[35] + x[18] <= 1
 clq714: x[54] + x[60] <= 1
 clq715: x[58] + x[80] <= 1
 clq716: x[84] + x[125] <= 1
 clq717: x[86] + x[102] <= 1
 clq718: x[2] + x[87] <= 1
 clq719: x[121] + x[106] <= 1
 clq720: x[122] + x[78] <= 1
 clq721: x[123] + x[71] <= 1
 clq722: x[141] + x[101] <= 1
 clq723: x[28] + x[23] <= 1
 clq724: x[36] + x[4] <= 1
 clq725: x[123] + x[53] + x[81] <= 1
 clq726: x[131] + x[62] <= 1
 clq727: x[65] + x[66] <= 1
 clq728: x[97] + x[105] <= 1
 clq729: x[128] + x[101] <= 1
 clq730: x[10] + x[45] <= 1
 clq731: x[143] + x[32] <= 1
 clq732: x[85] + x[44] <= 1
 clq733: x[92] + x[34] <= 1
 clq734: x[144] + x[150] <= 1
 clq735: x[113] + x[40] <= 1
 clq736: x[20] + x[119] <= 1
 clq737: x[39] + x[88] <= 1
 clq738: x[72] + x[32] + x[38] <= 1
 clq739: x[75] + x[11] <= 1
 clq740: x[33] + x[90] <= 1
 clq741: x[127] + x[100] <= 1
 clq742: x[126] + x[110] <= 1
 clq743: x[65] + x[30] <= 1
 clq744: x[37] + x[9] <= 1
 clq745: x[11] + x[89] <= 1
 clq746: x[139] + x[96] <= 1
 clq747: x[3] + x[122] <= 1
 clq748: x[5] + x[84] <= 1
 clq749: x[6] + x[90] <= 1
 clq750: x[136] + x[117] <= 1
 clq751: x[137] + x[51] <= 1
 clq752: x[14] + x[17] <= 1
 clq753: x[15] + x[43] <= 1
 clq754: x[21] + x[82] <= 1
 clq755: x[27] + x[26] <= 1
 clq756: x[31] + x[29] <= 1
 clq757: x[46] + x[125] <= 1
 clq758: x[49] + x[57] <= 1
 clq759: x[50] + x[111] <= 1
 clq760: x[67] + x[106] <= 1
 clq761: x[68] + x[129] <= 1
 clq762: x[76] + x[127] <= 1
 clq763: x[91] + x[62] <= 1
 clq764: x[93] + x[54] <= 1
 clq765: x[98] + x[59] <= 1
 clq766: x[103] + x[150] <= 1
 clq767: x[107] + x[25] <= 1
 clq768: x[109] + x[20] <= 1
 clq769: x[116] + x[19] <= 1
 clq770: x[134] + x[144] <= 1
 clq771: x[145] + x[61] <= 1
 clq772: x[24] + x[66] <= 1
 clq773: x[35] + x[73] <= 1
 clq774: x[48] + x[89] <= 1
 clq775: x[58] + x[85] <= 1
 clq776: x[69] + x[36] <= 1
 clq777: x[70] + x[28] <= 1
 clq778: x[77] + x[113] <= 1
 clq779: x[79] + x[121] <= 1
 clq780: x[86] + x[39] <= 1
 clq781: x[137] + x[87] <= 1
 clq782: x[141] + x[52] <= 1
 clq783: x[49] + x[148] <= 1
 clq784: x[60] + x[9] <= 1
 clq785: x[64] + x[135] <= 1
 clq786: x[80] + x[119] <= 1
 clq787: x[83] + x[118] <= 1
 clq788: x[50] + x[95] <= 1
 clq789: x[105] + x[22] <= 1
 clq790: x[97] + x[115] <= 1
 clq791: x[4] + x[71] <= 1
 clq792: x[99] + x[128] <= 1
 clq793: x[10] + x[132] <= 1
 clq794: x[143] + x[78] <= 1
 clq795: x[16] + x[146] <= 1
 clq796: x[18] + x[44] <= 1
 clq797: x[55] + x[110] <= 1
 clq798: x[65] + x[92] <= 1
 clq799: x[104] + x[149] <= 1
 clq800: x[114] + x[45] <= 1
 clq801: x[88] + x[7] <= 1
 clq802: x[8] + x[75] + x[40] <= 1
 clq803: x[66] + x[140] <= 1
 clq804: x[23] + x[110] <= 1
 clq805: x[41] + x[13] <= 1
 clq806: x[33] + x[126] <= 1
 clq807: x[67] + x[30] <= 1
 clq808: x[74] + x[63] <= 1
 clq809: x[19] + x[142] <= 1
 clq810: x[129] + x[42] <= 1
 clq811: x[35] + x[47] + x[132] <= 1
 clq812: x[125] + x[102] <= 1
 clq813: x[100] + x[34] <= 1
 clq814: x[25] + x[96] <= 1
 clq815: x[77] + x[120] <= 1
 clq816: x[130] + x[76] <= 1
 clq817: x[3] + x[84] <= 1
 clq818: x[5] + x[147] <= 1
 clq819: x[6] + x[114] <= 1
 clq820: x[136] + x[60] <= 1
 clq821: x[14] + x[89] <= 1
 clq822: x[15] + x[83] <= 1
 clq823: x[21] + x[92] <= 1
 clq824: x[27] + x[94] <= 1
 clq825: x[31] + x[90] <= 1
 clq826: x[46] + x[91] <= 1
 clq827: x[57] + x[28] <= 1
 clq828: x[68] + x[122] <= 1
 clq829: x[82] + x[85] <= 1
 clq830: x[93] + x[142] <= 1
 clq831: x[98] + x[95] <= 1
 clq832: x[103] + x[87] <= 1
 clq833: x[107] + x[56] <= 1
 clq834: x[109] + x[45] <= 1
 clq835: x[111] + x[18] <= 1
 clq836: x[112] + x[70] <= 1
 clq837: x[116] + x[126] <= 1
 clq838: x[124] + x[23] <= 1
 clq839: x[2] + x[121] <= 1
 clq840: x[131] + x[146] <= 1
 clq841: x[134] + x[128] <= 1
 clq842: x[138] + x[9] <= 1
 clq843: x[139] + x[145] <= 1
 clq844: x[24] + x[34] <= 1
 clq845: x[48] + x[61] <= 1
 clq846: x[116] + x[54] <= 1
 clq847: x[97] + x[58] <= 1
 clq848: x[145] + x[69] <= 1
 clq849: x[79] + x[26] <= 1
 clq850: x[86] + x[12] <= 1
 clq851: x[117] + x[53] <= 1
 clq852: x[141] + x[44] <= 1
 clq853: x[148] + x[20] <= 1
 clq854: x[29] + x[96] <= 1
 clq855: x[36] + x[106] <= 1
 clq856: x[37] + x[4] <= 1
 clq857: x[43] + x[105] <= 1
 clq858: x[27] + x[62] <= 1
 clq859: x[137] + x[91] + x[64] <= 1
 clq860: x[80] + x[150] <= 1
 clq861: x[115] + x[81] <= 1
 clq862: x[10] + x[140] <= 1
 clq863: x[143] + x[13] <= 1
 clq864: x[16] + x[44] <= 1
 clq865: x[70] + x[55] <= 1
 clq866: x[71] + x[39] <= 1
 clq867: x[144] + x[101] <= 1
 clq868: x[74] + x[104] <= 1
 clq869: x[113] + x[42] <= 1
 clq870: x[118] + x[63] <= 1
 clq871: x[5] + x[149] <= 1
 clq872: x[72] + x[135] <= 1
 clq873: x[88] + x[11] <= 1
 clq874: x[119] + x[7] <= 1
 clq875: x[133] + x[102] <= 1
 clq876: x[22] + x[73] <= 1
 clq877: x[33] + x[41] <= 1
 clq878: x[127] + x[73] <= 1
 clq879: x[30] + x[52] <= 1
 clq880: x[147] + x[59] <= 1
 clq881: x[9] + x[17] <= 1
 clq882: x[65] + x[51] <= 1
 clq883: x[133] + x[78] <= 1
 clq884: x[53] + x[120] <= 1
 clq885: x[1] + x[82] <= 1
 clq886: x[130] + x[69] <= 1
 clq887: x[3] + x[115] <= 1
 clq888: x[6] + x[149] <= 1
 clq889: x[136] + x[29] <= 1
 clq890: x[14] + x[113] <= 1
 clq891: x[15] + x[123] <= 1
 clq892: x[21] + x[143] <= 1
 clq893: x[31] + x[47] <= 1
 clq894: x[46] + x[124] <= 1
 clq895: x[49] + x[96] <= 1
 clq896: x[50] + x[51] <= 1
 clq897: x[57] + x[25] <= 1
 clq898: x[67] + x[83] <= 1
 clq899: x[68] + x[148] <= 1
 clq900: x[76] + x[37] <= 1
 clq901: x[93] + x[90] <= 1
 clq902: x[98] + x[108] <= 1
 clq903: x[99] + x[146] <= 1
 clq904: x[103] + x[100] <= 1
 clq905: x[107] + x[122] <= 1
 clq906: x[57] + x[109] <= 1
 clq907: x[1] + x[111] <= 1
 clq908: x[21] + x[112] <= 1
 clq909: x[2] + x[58] <= 1
 clq910: x[131] + x[141] <= 1
 clq911: x[134] + x[24] <= 1
 clq912: x[138] + x[12] <= 1
 clq913: x[139] + x[30] <= 1
 clq914: x[19] + x[34] <= 1
 clq915: x[48] + x[62] <= 1
 clq916: x[54] + x[94] <= 1
 clq917: x[61] + x[66] <= 1
 clq918: x[77] + x[122] <= 1
 clq919: x[79] + x[55] <= 1
 clq920: x[84] + x[60] <= 1
 clq921: x[86] + x[26] <= 1
 clq922: x[87] + x[17] <= 1
 clq923: x[117] + x[102] <= 1
 clq924: x[121] + x[85] <= 1
 clq925: x[125] + x[63] <= 1
 clq926: x[28] + x[59] <= 1
 clq927: x[36] + x[126] <= 1
 clq928: x[43] + x[127] <= 1
 clq929: x[80] + x[95] <= 1
 clq930: x[105] + x[135] <= 1
 clq931: x[4] + x[16] <= 1
 clq932: x[49] + x[128] <= 1
 clq933: x[129] + x[101] <= 1
 clq934: x[10] + x[81] <= 1
 clq935: x[18] + x[20] <= 1
 clq936: x[43] + x[71] <= 1
 clq937: x[92] + x[72] <= 1
 clq938: x[144] + x[56] <= 1
 clq939: x[10] + x[104] <= 1
 clq940: x[114] + x[39] <= 1
 clq941: x[118] + x[120] <= 1
 clq942: x[150] + x[142] <= 1
 clq943: x[129] + x[88] <= 1
 clq944: x[106] + x[140] <= 1
 clq945: x[37] + x[119] <= 1
 clq946: x[131] + x[8] <= 1
 clq947: x[85] + x[22] <= 1
 clq948: x[23] + x[45] <= 1
 clq949: x[76] + x[32] <= 1
 clq950: x[41] + x[75] <= 1
 clq951: x[145] + x[7] <= 1
 clq952: x[40] + x[78] <= 1
 clq953: x[42] + x[13] <= 1
 clq954: x[52] + x[11] <= 1
 clq955: x[110] + x[38] <= 1
 clq956: x[83] + x[89] <= 1
 clq957: x[130] + x[87] <= 1
 clq958: x[3] + x[148] <= 1
 clq959: x[5] + x[46] <= 1
 clq960: x[6] + x[111] <= 1
 clq961: x[136] + x[16] <= 1
 clq962: x[14] + x[55] <= 1
 clq963: x[15] + x[22] <= 1
 clq964: x[147] + x[95] <= 1
 clq965: x[27] + x[110] <= 1
 clq966: x[31] + x[24] <= 1
 clq967: x[33] + x[30] <= 1
 clq968: x[50] + x[109] <= 1
 clq969: x[67] + x[141] <= 1
 clq970: x[68] + x[71] <= 1
 clq971: x[74] + x[92] <= 1
 clq972: x[82] + x[127] <= 1
 clq973: x[93] + x[47] <= 1
 clq974: x[97] + x[38] <= 1
 clq975: x[98] + x[132] <= 1
 clq976: x[99] + x[117] <= 1
 clq977: x[103] + x[104] <= 1
 clq978: x[107] + x[28] <= 1
 clq979: x[108] + x[121] <= 1
 clq980: x[112] + x[64] <= 1
 clq981: x[116] + x[69] <= 1
 clq982: x[124] + x[150] <= 1
 clq983: x[2] + x[86] <= 1
 clq984: x[134] + x[123] <= 1
 clq985: x[138] + x[58] <= 1
 clq986: x[139] + x[128] <= 1
 clq987: x[19] + x[60] <= 1
 clq988: x[25] + x[41] <= 1
 clq989: x[35] + x[77] <= 1
 clq990: x[48] + x[144] <= 1
 clq991: x[54] + x[8] <= 1
 clq992: x[61] + x[18] <= 1
 clq993: x[70] + x[94] <= 1
 clq994: x[79] + x[53] <= 1
 clq995: x[84] + x[96] <= 1
 clq996: x[125] + x[62] <= 1
 clq997: x[29] + x[146] <= 1
 clq998: x[36] + x[88] <= 1
 clq999: x[65] + x[114] <= 1
 clq1000: x[80] + x[44] <= 1
 clq1001: x[105] + x[56] <= 1
 clq1002: x[115] + x[89] <= 1
 clq1003: x[4] + x[113] <= 1
 clq1004: x[143] + x[140] <= 1
 clq1005: x[26] + x[40] <= 1
 clq1006: x[81] + x[126] <= 1
 clq1007: x[118] + x[20] <= 1
 clq1008: x[149] + x[63] <= 1
 clq1009: x[39] + x[66] <= 1
 clq1010: x[72] + x[120] <= 1
 clq1011: x[106] + x[7] <= 1
 clq1012: x[119] + x[135] <= 1
 clq1013: x[133] + x[32] <= 1
 clq1014: x[23] + x[78] <= 1
 clq1015: x[75] + x[101] <= 1
 clq1016: x[90] + x[11] <= 1
 clq1017: x[12] + x[52] <= 1
 clq1018: x[59] + x[73] <= 1
 clq1019: x[9] + x[42] <= 1
 clq1020: x[142] + x[100] <= 1
 clq1021: x[17] + x[45] <= 1
 clq1022: x[102] + x[51] <= 1
 clq1023: x[13] + x[34] <= 1
Bounds
Binaries
 x[1] x[130] x[3] x[5] x[6] x[136] x[137] x[14] x[15] x[147] x[21] x[27]
 x[31] x[33] x[46] x[49] x[50] x[57] x[67] x[68] x[74] x[76] x[82] x[91]
 x[93] x[97] x[98] x[99] x[103] x[107] x[108] x[109] x[111] x[112] x[116]
 x[124] x[2] x[131] x[134] x[138] x[139] x[145] x[19] x[24] x[25] x[35]
 x[48] x[54] x[58] x[61] x[69] x[70] x[77] x[79] x[84] x[86] x[87] x[117]
 x[121] x[122] x[123] x[125] x[141] x[148] x[28] x[29] x[36] x[37] x[43]
 x[53] x[60] x[62] x[64] x[65] x[80] x[83] x[95] x[105] x[115] x[4] x[128]
 x[129] x[10] x[143] x[16] x[18] x[26] x[55] x[71] x[81] x[85] x[92] x[144]
 x[104] x[113] x[114] x[118] x[146] x[20] x[149] x[150] x[39] x[66] x[72]
 x[88] x[106] x[119] x[133] x[135] x[8] x[140] x[22] x[23] x[32] x[41]
 x[44] x[56] x[75] x[90] x[127] x[126] x[7] x[30] x[94] x[12] x[40] x[59]
 x[63] x[9] x[142] x[17] x[42] x[47] x[52] x[73] x[100] x[132] x[11] x[45]
 x[102] x[110] x[89] x[13] x[34] x[51] x[78] x[96] x[101] x[38] x[120]
End
