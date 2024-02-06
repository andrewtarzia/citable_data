%chk=cis_l2_la_FINAL.chk
%nprocshared=48
%mem=120GB
# SP PBE1PBE/def2svp empiricaldispersion=GD3BJ scrf=(solvent=dmso)

cis_l2_la

4 1
 Pd                -6.15340000   -3.29720000    0.15550000
 Pd                 6.15320000   -3.29750000   -0.15570000
 C                  0.10880000   -6.40500000    4.99640000
 C                  1.31530000   -5.98750000    4.44000000
 C                  1.32120000   -5.03680000    3.41100000
 C                  0.10460000   -4.52070000    2.95130000
 C                 -1.11320000   -4.96790000    3.47570000
 C                 -1.10000000   -5.91290000    4.51140000
 C                 -2.37710000   -4.49950000    2.87370000
 C                 -2.53000000   -3.18780000    2.40250000
 C                 -3.68440000   -2.82990000    1.72580000
 N                 -4.67540000   -3.70460000    1.49930000
 C                 -4.57570000   -4.94730000    1.99350000
 C                 -3.45590000   -5.37420000    2.68750000
 C                  2.57240000   -4.61340000    2.75200000
 C                  2.79170000   -3.27260000    2.40820000
 C                  3.91660000   -2.92800000    1.67980000
 N                  4.80870000   -3.84490000    1.27420000
 C                  4.64890000   -5.12380000    1.64110000
 C                  3.55750000   -5.54040000    2.38800000
 H                  0.11010000   -7.12900000    5.81370000
 H                  2.25820000   -6.38970000    4.81730000
 H                  0.10510000   -3.81650000    2.11640000
 H                 -2.03920000   -6.25750000    4.94920000
 H                 -1.75120000   -2.43860000    2.55050000
 H                 -3.81700000   -1.81900000    1.33480000
 H                 -5.40700000   -5.62680000    1.79580000
 H                 -3.41060000   -6.40920000    3.02820000
 H                  2.08780000   -2.49260000    2.70100000
 H                  4.10660000   -1.89160000    1.39330000
 H                  5.40090000   -5.83630000    1.29640000
 H                  3.45270000   -6.59900000    2.62910000
 C                 -0.10910000   -6.40550000   -4.99660000
 C                  1.09970000   -5.91340000   -4.51150000
 C                  1.11290000   -4.96840000   -3.47580000
 C                 -0.10490000   -4.52140000   -2.95130000
 C                 -1.32150000   -5.03740000   -3.41110000
 C                 -1.31560000   -5.98800000   -4.44020000
 C                 -2.57270000   -4.61380000   -2.75210000
 C                 -2.79170000   -3.27310000   -2.40830000
 C                 -3.91660000   -2.92820000   -1.67990000
 N                 -4.80880000   -3.84490000   -1.27430000
 C                 -4.64930000   -5.12390000   -1.64120000
 C                 -3.55800000   -5.54070000   -2.38820000
 C                  2.37670000   -4.50000000   -2.87370000
 C                  2.52950000   -3.18830000   -2.40230000
 C                  3.68380000   -2.83040000   -1.72560000
 N                  4.67510000   -3.70490000   -1.49940000
 C                  4.57540000   -4.94750000   -1.99380000
 C                  3.45570000   -5.37450000   -2.68780000
 H                 -0.11040000   -7.12940000   -5.81390000
 H                  2.03900000   -6.25790000   -4.94930000
 H                 -0.10550000   -3.81730000   -2.11640000
 H                 -2.25840000   -6.39030000   -4.81750000
 H                 -2.08770000   -2.49320000   -2.70110000
 H                 -4.10630000   -1.89180000   -1.39330000
 H                 -5.40150000   -5.83620000   -1.29660000
 H                 -3.45340000   -6.59930000   -2.62930000
 H                  1.75060000   -2.43920000   -2.55010000
 H                  3.81610000   -1.81950000   -1.33440000
 H                  5.40690000   -5.62700000   -1.79620000
 H                  3.41060000   -6.40940000   -3.02860000
 C                 -9.53180000   -3.76080000   -2.36800000
 C                 -8.42010000   -3.97880000   -1.60020000
 N                 -7.63000000   -2.94900000   -1.18760000
 C                 -7.89810000   -1.71230000   -1.56310000
 C                 -8.99960000   -1.38310000   -2.38130000
 C                 -9.28700000   -0.03480000   -2.76740000
 C                 -8.35340000    1.07090000   -2.46460000
 C                 -6.99710000    0.97750000   -2.81700000
 C                 -6.12420000    2.02670000   -2.57220000
 C                 -6.58320000    3.21420000   -1.97140000
 C                 -5.67190000    4.28830000   -1.76850000
 C                 -4.83860000    5.17230000   -1.66850000
 C                 -3.81670000    6.16290000   -1.61480000
 C                 -2.47260000    5.75780000   -1.65540000
 C                 -1.48010000    6.73250000   -1.65180000
 C                 -0.03170000    6.72400000   -1.68440000
 C                  0.94980000    5.73810000   -1.73480000
 C                  2.29820000    6.13060000   -1.74810000
 C                  3.32780000    5.14800000   -1.81660000
 C                  4.21200000    4.31340000   -1.90680000
 C                  5.24830000    3.35510000   -2.09560000
 C                  6.57940000    3.65240000   -1.74970000
 C                  7.59770000    2.74870000   -2.02230000
 C                  7.32600000    1.51880000   -2.64100000
 C                  8.43050000    0.61520000   -3.03190000
 C                  9.50480000    1.11290000   -3.74820000
 C                 10.54140000    0.27350000   -4.21050000
 C                 10.51710000   -1.08080000   -3.96990000
 C                  9.45920000   -1.63490000   -3.21240000
 C                  9.38980000   -3.01330000   -2.89080000
 C                  8.38360000   -3.48300000   -2.09080000
 N                  7.44020000   -2.64760000   -1.57580000
 C                  7.44880000   -1.36240000   -1.88280000
 C                  8.42150000   -0.78300000   -2.72410000
 C                  5.99220000    1.20960000   -2.95540000
 C                  4.96890000    2.10840000   -2.68650000
 C                  2.63910000    7.50410000   -1.71150000
 C                  1.66750000    8.49650000   -1.66390000
 C                  0.34380000    8.08030000   -1.65150000
 O                 -0.74440000    8.89740000   -1.60270000
 C                 -1.84130000    8.09180000   -1.60490000
 C                 -3.16040000    8.52180000   -1.55890000
 C                 -4.14340000    7.53940000   -1.56090000
 C                 -7.94110000    3.31020000   -1.61820000
 C                 -8.80890000    2.25340000   -1.86330000
 C                -10.46160000    0.20740000   -3.45640000
 C                -11.33570000   -0.84070000   -3.81830000
 C                -11.04640000   -2.14770000   -3.50000000
 C                 -9.87480000   -2.44420000   -2.76540000
 H                -10.16130000   -4.60270000   -2.66030000
 H                 -8.13180000   -4.97700000   -1.26710000
 H                 -7.23820000   -0.92830000   -1.18840000
 H                 -6.62630000    0.08150000   -3.32020000
 H                 -5.07600000    1.94490000   -2.86500000
 H                 -2.22320000    4.69630000   -1.68840000
 H                  0.68730000    4.67930000   -1.75960000
 H                  6.80690000    4.60620000   -1.27070000
 H                  8.62520000    2.99740000   -1.74760000
 H                  9.53240000    2.17650000   -3.99410000
 H                 11.36080000    0.71100000   -4.78480000
 H                 11.30650000   -1.73480000   -4.34520000
 H                 10.14860000   -3.70190000   -3.26570000
 H                  8.30670000   -4.53520000   -1.81350000
 H                  6.67750000   -0.73860000   -1.42640000
 H                  5.75510000    0.27020000   -3.46060000
 H                  3.94200000    1.86420000   -2.96430000
 H                  3.69360000    7.78350000   -1.71520000
 H                  1.93020000    9.55430000   -1.62720000
 H                 -3.40990000    9.58260000   -1.51570000
 H                 -5.19490000    7.82960000   -1.52700000
 H                 -8.31140000    4.22260000   -1.14750000
 H                 -9.85840000    2.33970000   -1.57310000
 H                -10.70020000    1.22920000   -3.75880000
 H                -12.24420000   -0.60400000   -4.37640000
 H                -11.71130000   -2.96090000   -3.79730000
 C                 -9.39020000   -3.01310000    2.89030000
 C                 -8.38400000   -3.48270000    2.09030000
 N                 -7.44020000   -2.64740000    1.57570000
 C                 -7.44860000   -1.36230000    1.88320000
 C                 -8.42130000   -0.78300000    2.72440000
 C                 -8.43020000    0.61520000    3.03250000
 C                 -7.32560000    1.51870000    2.64180000
 C                 -7.59720000    2.74860000    2.02310000
 C                 -6.57890000    3.65230000    1.75040000
 C                 -5.24780000    3.35480000    2.09620000
 C                 -4.21140000    4.31300000    1.90710000
 C                 -3.32730000    5.14760000    1.81660000
 C                 -2.29780000    6.13030000    1.74780000
 C                 -0.94930000    5.73790000    1.73470000
 C                  0.03210000    6.72390000    1.68410000
 C                  1.48050000    6.73250000    1.65170000
 C                  2.47300000    5.75780000    1.65560000
 C                  3.81710000    6.16300000    1.61510000
 C                  4.83900000    5.17250000    1.66890000
 C                  5.67230000    4.28830000    1.76900000
 C                  6.58330000    3.21400000    1.97180000
 C                  6.12420000    2.02670000    2.57290000
 C                  6.99690000    0.97730000    2.81760000
 C                  8.35310000    1.07040000    2.46500000
 C                  9.28670000   -0.03530000    2.76770000
 C                 10.46130000    0.20680000    3.45680000
 C                 11.33540000   -0.84120000    3.81850000
 C                 11.04620000   -2.14820000    3.50000000
 C                  9.87460000   -2.44470000    2.76540000
 C                  9.53180000   -3.76120000    2.36770000
 C                  8.42010000   -3.97920000    1.59980000
 N                  7.62990000   -2.94940000    1.18750000
 C                  7.89780000   -1.71280000    1.56320000
 C                  8.99930000   -1.38370000    2.38140000
 C                  8.80880000    2.25280000    1.86330000
 C                  7.94110000    3.30970000    1.61820000
 C                  4.14370000    7.53950000    1.56080000
 C                  3.16070000    8.52180000    1.55850000
 C                  1.84160000    8.09180000    1.60450000
 O                  0.74460000    8.89740000    1.60180000
 C                 -0.34350000    8.08010000    1.65070000
 C                 -1.66730000    8.49620000    1.66280000
 C                 -2.63880000    7.50380000    1.71070000
 C                 -4.96850000    2.10820000    2.68710000
 C                 -5.99180000    1.20940000    2.95610000
 C                 -9.50450000    1.11290000    3.74870000
 C                -10.54120000    0.27360000    4.21080000
 C                -10.51720000   -1.08070000    3.96970000
 C                 -9.45930000   -1.63480000    3.21230000
 H                -10.14910000   -3.70170000    3.26480000
 H                 -8.30720000   -4.53490000    1.81270000
 H                 -6.67690000   -0.73870000    1.42720000
 H                 -8.62470000    2.99740000    1.74840000
 H                 -6.80630000    4.60610000    1.27140000
 H                 -0.68680000    4.67920000    1.75990000
 H                  2.22370000    4.69630000    1.68880000
 H                  5.07600000    1.94520000    2.86600000
 H                  6.62610000    0.08150000    3.32110000
 H                 10.69980000    1.22860000    3.75920000
 H                 12.24390000   -0.60460000    4.37660000
 H                 11.71110000   -2.96140000    3.79710000
 H                 10.16140000   -4.60320000    2.65970000
 H                  8.13200000   -4.97740000    1.26640000
 H                  7.23770000   -0.92880000    1.18880000
 H                  9.85830000    2.33880000    1.57290000
 H                  8.31150000    4.22190000    1.14730000
 H                  5.19520000    7.82970000    1.52690000
 H                  3.41010000    9.58270000    1.51500000
 H                 -1.93000000    9.55400000    1.62580000
 H                 -3.69330000    7.78310000    1.71410000
 H                 -3.94150000    1.86380000    2.96480000
 H                 -5.75470000    0.27000000    3.46120000
 H                 -9.53190000    2.17640000    3.99490000
 H                -11.36070000    0.71100000    4.78510000
 H                -11.30680000   -1.73460000    4.34480000

 1 12 1.0 42 1.0 65 1.0 140 1.0
 2 18 1.0 48 1.0 94 1.0 169 1.0
 3 4 1.5 8 1.5 21 1.0
 4 5 1.5 22 1.0
 5 6 1.5 15 1.0
 6 7 1.5 23 1.0
 7 8 1.5 9 1.0
 8 24 1.0
 9 10 1.5 14 1.5
 10 11 2.0 25 1.0
 11 12 1.5 26 1.0
 12 13 1.5
 13 14 2.0 27 1.0
 14 28 1.0
 15 16 1.5 20 1.5
 16 17 2.0 29 1.0
 17 18 1.5 30 1.0
 18 19 1.5
 19 20 1.5 31 1.0
 20 32 1.0
 21
 22
 23
 24
 25
 26
 27
 28
 29
 30
 31
 32
 33 34 1.5 38 1.5 51 1.0
 34 35 1.5 52 1.0
 35 36 1.5 45 1.0
 36 37 1.5 53 1.0
 37 38 1.5 39 1.0
 38 54 1.0
 39 40 1.5 44 1.5
 40 41 2.0 55 1.0
 41 42 1.5 56 1.0
 42 43 1.5
 43 44 1.5 57 1.0
 44 58 1.0
 45 46 1.5 50 1.5
 46 47 2.0 59 1.0
 47 48 1.5 60 1.0
 48 49 1.5
 49 50 2.0 61 1.0
 50 62 1.0
 51
 52
 53
 54
 55
 56
 57
 58
 59
 60
 61
 62
 63 64 2.0 111 1.5 112 1.0
 64 65 1.5 113 1.0
 65 66 2.0
 66 67 1.5 114 1.0
 67 68 1.5 111 1.5
 68 69 1.0 108 2.0
 69 70 1.5 107 1.5
 70 71 1.5 115 1.0
 71 72 1.5 116 1.0
 72 73 1.5 106 1.5
 73 74 3.0
 74 75 1.5
 75 76 1.5 105 1.5
 76 77 1.5 117 1.0
 77 78 1.0 103 1.5
 78 79 1.5 101 1.5
 79 80 1.5 118 1.0
 80 81 1.5 99 1.5
 81 82 3.0
 82 83 1.5
 83 84 1.5 98 1.5
 84 85 1.5 119 1.0
 85 86 1.5 120 1.0
 86 87 1.0 97 1.5
 87 88 2.0 96 1.5
 88 89 1.5 121 1.0
 89 90 2.0 122 1.0
 90 91 1.5 123 1.0
 91 92 1.5 96 1.5
 92 93 2.0 124 1.0
 93 94 1.5 125 1.0
 94 95 2.0
 95 96 1.5 126 1.0
 96
 97 98 1.5 127 1.0
 98 128 1.0
 99 100 1.5 129 1.0
 100 101 1.5 130 1.0
 101 102 1.0
 102 103 1.0
 103 104 1.5
 104 105 1.5 131 1.0
 105 132 1.0
 106 107 1.5 133 1.0
 107 134 1.0
 108 109 1.5 135 1.0
 109 110 2.0 136 1.0
 110 111 1.5 137 1.0
 111
 112
 113
 114
 115
 116
 117
 118
 119
 120
 121
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
 138 139 2.0 186 1.5 187 1.0
 139 140 1.5 188 1.0
 140 141 2.0
 141 142 1.5 189 1.0
 142 143 1.5 186 1.5
 143 144 1.0 183 2.0
 144 145 1.5 182 1.5
 145 146 1.5 190 1.0
 146 147 1.5 191 1.0
 147 148 1.5 181 1.5
 148 149 3.0
 149 150 1.5
 150 151 1.5 180 1.5
 151 152 1.5 192 1.0
 152 153 1.0 178 1.5
 153 154 1.5 176 1.5
 154 155 1.5 193 1.0
 155 156 1.5 174 1.5
 156 157 3.0
 157 158 1.5
 158 159 1.5 173 1.5
 159 160 1.5 194 1.0
 160 161 1.5 195 1.0
 161 162 1.0 172 1.5
 162 163 2.0 171 1.5
 163 164 1.5 196 1.0
 164 165 2.0 197 1.0
 165 166 1.5 198 1.0
 166 167 1.5 171 1.5
 167 168 2.0 199 1.0
 168 169 1.5 200 1.0
 169 170 2.0
 170 171 1.5 201 1.0
 171
 172 173 1.5 202 1.0
 173 203 1.0
 174 175 1.5 204 1.0
 175 176 1.5 205 1.0
 176 177 1.0
 177 178 1.0
 178 179 1.5
 179 180 1.5 206 1.0
 180 207 1.0
 181 182 1.5 208 1.0
 182 209 1.0
 183 184 1.5 210 1.0
 184 185 2.0 211 1.0
 185 186 1.5 212 1.0
 186
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
 211
 212

