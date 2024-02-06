%chk=trans_l1_lb_FINAL.chk
%nprocshared=48
%mem=120GB
# SP PBE1PBE/def2svp empiricaldispersion=GD3BJ scrf=(solvent=dmso)

trans_l1_lb

4 1
 Pd                 6.14320000    1.94780000    1.33460000
 Pd                -6.69200000    1.55030000   -0.27530000
 C                  4.66490000    0.11660000   -0.17240000
 N                  5.68160000    0.00480000    0.70460000
 C                  5.95650000   -1.20950000    1.20600000
 C                  5.18050000   -2.32790000    0.90980000
 C                  4.06460000   -2.20680000    0.06610000
 C                  3.09760000   -3.28730000   -0.21420000
 C                  3.49580000   -4.62440000   -0.42440000
 C                  2.57920000   -5.61800000   -0.76290000
 C                  1.24670000   -5.24140000   -0.86460000
 O                  0.19970000   -6.05460000   -1.18410000
 C                 -0.92760000   -5.29240000   -1.13080000
 C                 -2.22220000   -5.71760000   -1.39960000
 C                 -3.24230000   -4.78060000   -1.24970000
 C                 -2.98310000   -3.45420000   -0.84030000
 C                 -4.05990000   -2.46840000   -0.62700000
 C                 -3.88010000   -1.44600000    0.31840000
 C                 -4.80870000   -0.43480000    0.43430000
 N                 -5.92420000   -0.37870000   -0.32060000
 C                 -6.16380000   -1.38870000   -1.17050000
 C                 -5.27070000   -2.44290000   -1.33840000
 C                 -1.65890000   -3.04810000   -0.61570000
 C                 -0.62870000   -3.96880000   -0.76410000
 C                  0.81020000   -3.92870000   -0.61400000
 C                  1.74080000   -2.94390000   -0.30160000
 C                  3.85520000   -0.94610000   -0.51360000
 H                  4.47740000    1.10450000   -0.59920000
 H                  6.80010000   -1.27940000    1.89560000
 H                  5.42590000   -3.27680000    1.38980000
 H                  4.55460000   -4.88330000   -0.35910000
 H                  2.89430000   -6.64520000   -0.94690000
 H                 -2.43240000   -6.74610000   -1.69310000
 H                 -4.27300000   -5.09950000   -1.41720000
 H                 -3.00790000   -1.42390000    0.97090000
 H                 -4.64120000    0.37080000    1.15130000
 H                 -7.08290000   -1.33360000   -1.75680000
 H                 -5.51650000   -3.21630000   -2.06700000
 H                 -1.42920000   -2.01650000   -0.34560000
 H                  1.41210000   -1.92200000   -0.10360000
 H                  3.05100000   -0.78280000   -1.23100000
 C                  4.21460000    3.98800000    2.08680000
 N                  5.39510000    3.92150000    1.44060000
 C                  5.68870000    4.91140000    0.58290000
 C                  4.81690000    5.95530000    0.31090000
 C                  3.55740000    5.99870000    0.92500000
 C                  2.53200000    6.99780000    0.57200000
 C                  2.88610000    8.31320000    0.19980000
 C                  1.93190000    9.25590000   -0.17130000
 C                  0.60640000    8.84080000   -0.17690000
 O                 -0.46600000    9.59550000   -0.54090000
 C                 -1.56340000    8.79720000   -0.43660000
 C                 -2.86780000    9.15850000   -0.74880000
 C                 -3.84480000    8.17680000   -0.61010000
 C                 -3.53630000    6.87720000   -0.15280000
 C                 -4.56660000    5.82700000   -0.05470000
 C                 -5.61200000    5.70160000   -0.98020000
 C                 -6.43630000    4.58560000   -0.94840000
 N                 -6.29140000    3.60600000   -0.04450000
 C                 -5.36580000    3.77270000    0.92130000
 C                 -4.50160000    4.84930000    0.94950000
 C                 -2.21280000    6.54730000    0.16940000
 C                 -1.21970000    7.51020000    0.01530000
 C                  0.21780000    7.53840000    0.18580000
 C                  1.18430000    6.61330000    0.57020000
 C                  3.29250000    4.99240000    1.86490000
 H                  3.98770000    3.19030000    2.79240000
 H                  6.64370000    4.85960000    0.06110000
 H                  5.11290000    6.69830000   -0.43090000
 H                  3.93630000    8.61120000    0.22050000
 H                  2.21040000   10.27300000   -0.44940000
 H                 -3.11230000   10.16420000   -1.09260000
 H                 -4.88080000    8.43020000   -0.84440000
 H                 -5.75460000    6.43340000   -1.77650000
 H                 -7.21350000    4.46070000   -1.70210000
 H                 -5.30220000    2.99890000    1.68790000
 H                 -3.76640000    4.90950000    1.75230000
 H                 -1.95390000    5.53830000    0.49550000
 H                  0.89520000    5.59220000    0.82500000
 H                  2.35930000    4.97470000    2.42870000
 C                  6.35690000    2.26350000    5.48330000
 C                  6.25130000    2.49830000    4.14160000
 N                  5.66110000    1.60050000    3.29410000
 C                  5.14450000    0.48690000    3.78040000
 C                  5.22120000    0.13930000    5.14700000
 C                  4.75590000   -1.12050000    5.63790000
 C                  4.07210000   -2.10460000    4.77140000
 C                  2.90200000   -1.76870000    4.07290000
 C                  2.19460000   -2.73660000    3.37540000
 C                  2.64200000   -4.07080000    3.35560000
 C                  1.80760000   -5.06560000    2.77500000
 C                  0.97910000   -5.85430000    2.35930000
 C                 -0.13550000   -6.64680000    1.93760000
 N                 -1.32560000   -6.06740000    2.12910000
 C                 -2.43060000   -6.72510000    1.76530000
 C                 -3.64760000   -6.01560000    2.02840000
 C                 -4.56550000   -5.27770000    2.33590000
 C                 -5.52120000   -4.32820000    2.79940000
 C                 -5.36730000   -3.80270000    4.09540000
 C                 -6.21560000   -2.80680000    4.55600000
 C                 -7.23340000   -2.30050000    3.73680000
 C                 -8.07020000   -1.18110000    4.22430000
 C                 -8.78770000   -1.28770000    5.40020000
 C                 -9.58040000   -0.22060000    5.87870000
 C                 -9.66200000    0.96840000    5.19080000
 C                 -8.92760000    1.13220000    3.99350000
 C                 -8.94180000    2.32590000    3.23150000
 C                 -8.19740000    2.42400000    2.08800000
 N                 -7.40340000    1.40050000    1.64920000
 C                 -7.36770000    0.27440000    2.34000000
 C                 -8.12150000    0.05880000    3.51180000
 C                 -7.41080000   -2.85290000    2.45840000
 C                 -6.56930000   -3.85570000    1.99220000
 C                 -2.39450000   -8.00120000    1.18290000
 C                 -1.15240000   -8.59870000    0.98840000
 C                  0.00180000   -7.92350000    1.37060000
 C                  3.84700000   -4.39600000    3.99970000
 C                  4.54560000   -3.42290000    4.70270000
 C                  4.97110000   -1.42230000    6.97070000
 C                  5.59510000   -0.50840000    7.84640000
 C                  6.03000000    0.71620000    7.39540000
 C                  5.86770000    1.05370000    6.03270000
 H                  6.84560000    3.00090000    6.12190000
 H                  6.64820000    3.40700000    3.68820000
 H                  4.67240000   -0.19510000    3.07460000
 H                  2.51110000   -0.74970000    4.12130000
 H                  1.26340000   -2.47850000    2.86840000
 H                 -4.55840000   -4.17040000    4.72920000
 H                 -6.06830000   -2.38740000    5.55360000
 H                 -8.75680000   -2.22450000    5.96070000
 H                -10.14180000   -0.35260000    6.80620000
 H                -10.28170000    1.78950000    5.55620000
 H                 -9.55070000    3.17400000    3.54850000
 H                 -8.20880000    3.32480000    1.47640000
 H                 -6.72670000   -0.52160000    1.96840000
 H                 -8.22010000   -2.49310000    1.81850000
 H                 -6.70900000   -4.26660000    0.99070000
 H                 -3.31950000   -8.50200000    0.89470000
 H                 -1.08490000   -9.59070000    0.53730000
 H                  0.99160000   -8.35850000    1.22730000
 H                  4.21310000   -5.42380000    3.97390000
 H                  5.46930000   -3.68910000    5.22110000
 H                  4.62680000   -2.38300000    7.35900000
 H                  5.72700000   -0.78310000    8.89510000
 H                  6.51410000    1.42650000    8.06830000
 C                  9.39650000    3.45470000   -0.81650000
 C                  8.35170000    3.21070000    0.03010000
 N                  7.38760000    2.28340000   -0.25230000
 C                  7.50960000    1.53010000   -1.32880000
 C                  8.60190000    1.63800000   -2.21910000
 C                  8.76500000    0.75410000   -3.33180000
 C                  7.90550000   -0.43760000   -3.48430000
 C                  7.77050000   -1.35360000   -2.42820000
 C                  6.95130000   -2.46690000   -2.55160000
 C                  6.23970000   -2.68830000   -3.74240000
 C                  5.29170000   -3.74600000   -3.84730000
 C                  4.38770000   -4.55480000   -3.94430000
 C                  3.16110000   -5.28930000   -4.04830000
 N                  2.07590000   -4.51480000   -4.12650000
 C                  0.87330000   -5.09380000   -4.20170000
 C                 -0.21930000   -4.17360000   -4.26970000
 C                 -1.03570000   -3.27330000   -4.33130000
 C                 -1.91350000   -2.16420000   -4.47250000
 C                 -3.29490000   -2.28380000   -4.25070000
 C                 -4.13160000   -1.19860000   -4.46880000
 C                 -3.61820000    0.02510000   -4.92820000
 C                 -4.55360000    1.09180000   -5.33960000
 C                 -4.51760000    1.56620000   -6.63890000
 C                 -5.51280000    2.43390000   -7.13600000
 C                 -6.57900000    2.80610000   -6.34950000
 C                 -6.63380000    2.37930000   -5.00360000
 C                 -7.69980000    2.72020000   -4.13600000
 C                 -7.64590000    2.35730000   -2.82060000
 N                 -6.55970000    1.71920000   -2.29140000
 C                 -5.58300000    1.31380000   -3.08150000
 C                 -5.58660000    1.56610000   -4.47400000
 C                 -2.23560000    0.15180000   -5.11950000
 C                 -1.39110000   -0.92380000   -4.88600000
 C                  0.70370000   -6.48750000   -4.21390000
 C                  1.83930000   -7.28650000   -4.13860000
 C                  3.09440000   -6.69050000   -4.05210000
 C                  6.40400000   -1.79400000   -4.81530000
 C                  7.23040000   -0.68650000   -4.68630000
 C                  9.77960000    1.01750000   -4.23360000
 C                 10.67380000    2.09290000   -4.04240000
 C                 10.58230000    2.90110000   -2.93180000
 C                  9.54490000    2.68500000   -1.99650000
 H                 10.12070000    4.23230000   -0.56900000
 H                  8.24920000    3.75100000    0.97290000
 H                  6.72600000    0.79920000   -1.52400000
 H                  8.32330000   -1.19580000   -1.49900000
 H                  6.84370000   -3.16500000   -1.71950000
 H                 -3.70440000   -3.24450000   -3.93540000
 H                 -5.20970000   -1.31570000   -4.33780000
 H                 -3.73270000    1.21340000   -7.31130000
 H                 -5.45290000    2.77340000   -8.17220000
 H                 -7.38190000    3.43080000   -6.74530000
 H                 -8.56110000    3.27570000   -4.50970000
 H                 -8.45990000    2.58790000   -2.13080000
 H                 -4.73950000    0.79900000   -2.61640000
 H                 -1.82240000    1.10240000   -5.46380000
 H                 -0.31600000   -0.82280000   -5.04440000
 H                 -0.29570000   -6.92000000   -4.27160000
 H                  1.74690000   -8.37450000   -4.14230000
 H                  4.00540000   -7.28630000   -3.98510000
 H                  5.85620000   -1.96410000   -5.74380000
 H                  7.33060000    0.01520000   -5.51730000
 H                  9.91750000    0.35270000   -5.08900000
 H                 11.46790000    2.25870000   -4.77340000
 H                 11.29900000    3.70620000   -2.75970000

 1 4 1.0 83 1.0 148 1.0
 2 20 1.0 109 1.0 174 1.0
 3 4 1.5 27 2.0 28 1.0
 4 5 1.5
 5 6 1.5 29 1.0
 6 7 1.5 30 1.0
 7 8 1.0 27 1.5
 8 9 1.5 26 1.5
 9 10 1.5 31 1.0
 10 11 1.5 32 1.0
 11 12 1.0 25 1.5
 12 13 1.0
 13 14 1.5 24 1.5
 14 15 1.5 33 1.0
 15 16 1.5 34 1.0
 16 17 1.0 23 1.5
 17 18 1.5 22 1.5
 18 19 2.0 35 1.0
 19 20 1.5 36 1.0
 20 21 1.5
 21 22 1.5 37 1.0
 22 38 1.0
 23 24 1.5 39 1.0
 24 25 1.5
 25 26 1.5
 26 40 1.0
 27 41 1.0
 28
 29
 30
 31
 32
 33
 34
 35
 36
 37
 38
 39
 40
 41
 42 43 1.5 66 2.0 67 1.0
 43 44 1.5
 44 45 1.5 68 1.0
 45 46 1.5 69 1.0
 46 47 1.0 66 1.5
 47 48 1.5 65 1.5
 48 49 1.5 70 1.0
 49 50 1.5 71 1.0
 50 51 1.0 64 1.5
 51 52 1.0
 52 53 1.5 63 1.5
 53 54 1.5 72 1.0
 54 55 1.5 73 1.0
 55 56 1.0 62 1.5
 56 57 1.5 61 1.5
 57 58 1.5 74 1.0
 58 59 1.5 75 1.0
 59 60 1.5
 60 61 2.0 76 1.0
 61 77 1.0
 62 63 1.5 78 1.0
 63 64 1.0
 64 65 1.5
 65 79 1.0
 66 80 1.0
 67
 68
 69
 70
 71
 72
 73
 74
 75
 76
 77
 78
 79
 80
 81 82 2.0 122 1.5 123 1.0
 82 83 1.5 124 1.0
 83 84 2.0
 84 85 1.5 125 1.0
 85 86 1.5 122 1.5
 86 87 1.0 119 2.0
 87 88 1.5 118 1.5
 88 89 1.5 126 1.0
 89 90 1.5 127 1.0
 90 91 1.5 117 1.5
 91 92 3.0
 92 93 1.5
 93 94 1.5 116 1.5
 94 95 1.5
 95 96 1.5 114 1.5
 96 97 3.0
 97 98 1.5
 98 99 1.5 113 1.5
 99 100 1.5 128 1.0
 100 101 1.5 129 1.0
 101 102 1.0 112 1.5
 102 103 2.0 111 1.5
 103 104 1.5 130 1.0
 104 105 2.0 131 1.0
 105 106 1.5 132 1.0
 106 107 1.5 111 1.5
 107 108 2.0 133 1.0
 108 109 1.5 134 1.0
 109 110 2.0
 110 111 1.5 135 1.0
 111
 112 113 1.5 136 1.0
 113 137 1.0
 114 115 1.5 138 1.0
 115 116 1.5 139 1.0
 116 140 1.0
 117 118 1.5 141 1.0
 118 142 1.0
 119 120 1.5 143 1.0
 120 121 2.0 144 1.0
 121 122 1.5 145 1.0
 122
 123
 124
 125
 126
 127
 128
 129
 130
 131
 132
 133
 134
 135
 136
 137
 138
 139
 140
 141
 142
 143
 144
 145
 146 147 2.0 187 1.5 188 1.0
 147 148 1.5 189 1.0
 148 149 2.0
 149 150 1.5 190 1.0
 150 151 1.5 187 1.5
 151 152 1.0 184 2.0
 152 153 1.5 183 1.5
 153 154 1.5 191 1.0
 154 155 1.5 192 1.0
 155 156 1.5 182 1.5
 156 157 3.0
 157 158 1.5
 158 159 1.5 181 1.5
 159 160 1.5
 160 161 1.5 179 1.5
 161 162 3.0
 162 163 1.5
 163 164 1.5 178 1.5
 164 165 1.5 193 1.0
 165 166 1.5 194 1.0
 166 167 1.0 177 1.5
 167 168 2.0 176 1.5
 168 169 1.5 195 1.0
 169 170 2.0 196 1.0
 170 171 1.5 197 1.0
 171 172 1.5 176 1.5
 172 173 2.0 198 1.0
 173 174 1.5 199 1.0
 174 175 2.0
 175 176 1.5 200 1.0
 176
 177 178 1.5 201 1.0
 178 202 1.0
 179 180 1.5 203 1.0
 180 181 1.5 204 1.0
 181 205 1.0
 182 183 1.5 206 1.0
 183 207 1.0
 184 185 1.5 208 1.0
 185 186 2.0 209 1.0
 186 187 1.5 210 1.0
 187
 188
 189
 190
 191
 192
 193
 194
 195
 196
 197
 198
 199
 200
 201
 202
 203
 204
 205
 206
 207
 208
 209
 210

