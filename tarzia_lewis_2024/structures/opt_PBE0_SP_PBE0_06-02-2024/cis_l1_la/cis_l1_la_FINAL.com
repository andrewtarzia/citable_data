%chk=cis_l1_la_FINAL.chk
%nprocshared=48
%mem=120GB
# SP PBE1PBE/def2svp empiricaldispersion=GD3BJ scrf=(solvent=dmso)

cis_l1_la

4 1
 Pd                -7.34220000   -2.73770000    0.26500000
 Pd                 7.34170000   -2.73830000   -0.26410000
 C                 -6.14110000   -4.45300000    2.31580000
 N                 -5.97880000   -3.29240000    1.66370000
 C                 -4.86210000   -2.57850000    1.86950000
 C                 -3.85400000   -3.01840000    2.71020000
 C                 -3.98200000   -4.24750000    3.37530000
 C                 -2.89830000   -4.79900000    4.20950000
 C                 -3.19830000   -5.60750000    5.32830000
 C                 -2.20380000   -6.20490000    6.09460000
 C                 -0.89170000   -5.98260000    5.69970000
 O                  0.21370000   -6.51840000    6.28210000
 C                  1.28540000   -6.08150000    5.56680000
 C                  2.60890000   -6.43110000    5.79650000
 C                  3.55720000   -5.90220000    4.92690000
 C                  3.19920000   -5.03070000    3.87610000
 C                  4.22850000   -4.48990000    2.96880000
 C                  4.18430000   -3.15620000    2.53770000
 C                  5.13510000   -2.68880000    1.64800000
 N                  6.10350000   -3.48080000    1.16270000
 C                  6.17890000   -4.75360000    1.57520000
 C                  5.27370000   -5.28560000    2.47950000
 C                  1.85720000   -4.68700000    3.67390000
 C                  0.89110000   -5.23050000    4.51770000
 C                 -0.55360000   -5.17110000    4.60080000
 C                 -1.56320000   -4.56720000    3.85470000
 C                 -5.17680000   -4.95410000    3.17300000
 H                 -7.06270000   -5.00570000    2.12540000
 H                 -4.77510000   -1.62900000    1.33790000
 H                 -2.97190000   -2.39160000    2.84450000
 H                 -4.24080000   -5.76060000    5.61380000
 H                 -2.44080000   -6.82640000    6.95870000
 H                  2.88840000   -7.09230000    6.61750000
 H                  4.60970000   -6.15030000    5.07950000
 H                  3.42000000   -2.47060000    2.90540000
 H                  5.12680000   -1.65040000    1.30970000
 H                  6.97780000   -5.36450000    1.15010000
 H                  5.36530000   -6.33420000    2.76540000
 H                  1.56810000   -4.03640000    2.84600000
 H                 -1.31450000   -3.96340000    2.98020000
 H                 -5.34870000   -5.92170000    3.64540000
 C                 -5.13340000   -2.69050000   -1.64460000
 N                 -6.10330000   -3.48140000   -1.16060000
 C                 -6.18000000   -4.75400000   -1.57370000
 C                 -5.27480000   -5.28660000   -2.47770000
 C                 -4.22840000   -4.49200000   -2.96610000
 C                 -3.19940000   -5.03340000   -3.87330000
 C                 -3.55760000   -5.90490000   -4.92400000
 C                 -2.60950000   -6.43410000   -5.79360000
 C                 -1.28590000   -6.08490000   -5.56390000
 O                 -0.21430000   -6.52200000   -6.27930000
 C                  0.89120000   -5.98620000   -5.69710000
 C                  2.20330000   -6.20850000   -6.09210000
 C                  3.19780000   -5.61090000   -5.32600000
 C                  2.89800000   -4.80220000   -4.20730000
 C                  3.98170000   -4.25020000   -3.37340000
 C                  3.85320000   -3.02110000   -2.70850000
 C                  4.86120000   -2.58070000   -1.86780000
 N                  5.97810000   -3.29410000   -1.66210000
 C                  6.14090000   -4.45470000   -2.31410000
 C                  5.17670000   -4.95630000   -3.17120000
 C                  1.56290000   -4.57050000   -3.85240000
 C                  0.55330000   -5.17460000   -4.59830000
 C                 -0.89140000   -5.23380000   -4.51500000
 C                 -1.85720000   -4.69000000   -3.67120000
 C                 -4.18280000   -3.15860000   -2.53420000
 H                 -5.12390000   -1.65230000   -1.30570000
 H                 -6.97980000   -5.36410000   -1.14950000
 H                 -5.36750000   -6.33490000   -2.76410000
 H                 -4.61010000   -6.15270000   -5.07660000
 H                 -2.88930000   -7.09540000   -6.61440000
 H                  2.44030000   -6.83020000   -6.95610000
 H                  4.24030000   -5.76390000   -5.61160000
 H                  2.97080000   -2.39470000   -2.84270000
 H                  4.77390000   -1.63130000   -1.33630000
 H                  7.06280000   -5.00700000   -2.12370000
 H                  5.34900000   -5.92380000   -3.64350000
 H                  1.31430000   -3.96650000   -2.97800000
 H                 -1.56780000   -4.03950000   -2.84340000
 H                 -3.41750000   -2.47370000   -2.90100000
 C                -10.68710000   -2.78440000   -2.35390000
 C                 -9.65800000   -3.13500000   -1.52260000
 N                 -8.72800000   -2.22400000   -1.12220000
 C                 -8.78100000   -0.98070000   -1.56200000
 C                 -9.78120000   -0.52450000   -2.44650000
 C                 -9.83330000    0.83010000   -2.90560000
 C                 -8.76050000    1.79760000   -2.59040000
 C                 -7.41550000    1.49800000   -2.86180000
 C                 -6.41440000    2.42500000   -2.60980000
 C                 -6.72940000    3.69180000   -2.08400000
 C                 -5.70260000    4.65610000   -1.88370000
 C                 -4.80730000    5.47640000   -1.77650000
 C                 -3.76030000    6.43910000   -1.70690000
 C                 -2.41980000    6.02070000   -1.71540000
 C                 -1.41790000    6.98720000   -1.68760000
 C                  0.03260000    6.97090000   -1.68400000
 C                  1.01470000    5.98260000   -1.71210000
 C                  2.36340000    6.37400000   -1.69560000
 C                  3.41230000    5.41020000   -1.74970000
 C                  4.35360000    4.63960000   -1.83500000
 C                  5.49400000    3.80510000   -2.01570000
 C                  6.76890000    4.25770000   -1.62700000
 C                  7.90090000    3.50500000   -1.90730000
 C                  7.80370000    2.27490000   -2.57590000
 C                  9.02850000    1.55420000   -2.99060000
 C                 10.01680000    2.23190000   -3.68280000
 C                 11.17040000    1.57560000   -4.16360000
 C                 11.35330000    0.22700000   -3.96360000
 C                 10.38980000   -0.50560000   -3.23210000
 C                 10.53300000   -1.88710000   -2.95030000
 C                  9.60680000   -2.53010000   -2.17440000
 N                  8.54000000   -1.86730000   -1.64980000
 C                  8.35200000   -0.58760000   -1.91930000
 C                  9.23080000    0.16140000   -2.72880000
 C                  6.52670000    1.80270000   -2.92210000
 C                  5.38950000    2.55170000   -2.64710000
 C                  2.70550000    7.74750000   -1.64810000
 C                  1.73600000    8.74110000   -1.62170000
 C                  0.41240000    8.32610000   -1.64140000
 O                 -0.67140000    9.14870000   -1.61960000
 C                 -1.77120000    8.34940000   -1.65010000
 C                 -3.08650000    8.79220000   -1.63850000
 C                 -4.07650000    7.81890000   -1.66360000
 C                 -8.07550000    3.99120000   -1.80530000
 C                 -9.07150000    3.05720000   -2.05690000
 C                -10.91840000    1.21510000   -3.67210000
 C                -11.92890000    0.29900000   -4.03660000
 C                -11.86830000   -1.01770000   -3.64190000
 C                -10.79560000   -1.45390000   -2.82990000
 H                -11.43130000   -3.53020000   -2.63730000
 H                 -9.55200000   -4.14710000   -1.12970000
 H                 -8.01830000   -0.29390000   -1.19160000
 H                 -7.15140000    0.53680000   -3.30920000
 H                 -5.37550000    2.18390000   -2.84220000
 H                 -2.17950000    4.95690000   -1.74220000
 H                  0.75140000    4.92410000   -1.74350000
 H                  6.86140000    5.21640000   -1.11400000
 H                  8.88240000    3.87670000   -1.60480000
 H                  9.88350000    3.29510000   -3.89330000
 H                 11.91520000    2.15100000   -4.71760000
 H                 12.23540000   -0.28640000   -4.35090000
 H                 11.39190000   -2.43840000   -3.33570000
 H                  9.69440000   -3.58850000   -1.92560000
 H                  7.48710000   -0.10510000   -1.45960000
 H                  6.42340000    0.85870000   -3.46260000
 H                  4.40880000    2.18500000   -2.95560000
 H                  3.76050000    8.02400000   -1.62780000
 H                  1.99850000    9.79870000   -1.57760000
 H                 -3.32630000    9.85550000   -1.60250000
 H                 -5.12630000    8.11700000   -1.65530000
 H                 -8.33450000    4.96740000   -1.39130000
 H                -10.11100000    3.30260000   -1.82780000
 H                -10.97730000    2.24450000   -4.03190000
 H                -12.75920000    0.64460000   -4.65610000
 H                -12.63970000   -1.73100000   -3.93830000
 C                -10.53460000   -1.88470000    2.94930000
 C                 -9.60820000   -2.52830000    2.17400000
 N                 -8.54110000   -1.86590000    1.64960000
 C                 -8.35300000   -0.58610000    1.91860000
 C                 -9.23200000    0.16340000    2.72740000
 C                 -9.02960000    1.55630000    2.98870000
 C                 -7.80440000    2.27650000    2.57430000
 C                 -7.90100000    3.50630000    1.90510000
 C                 -6.76870000    4.25870000    1.62500000
 C                 -5.49410000    3.80590000    2.01460000
 C                 -4.35340000    4.64020000    1.83410000
 C                 -3.41210000    5.41070000    1.74900000
 C                 -2.36320000    6.37440000    1.69500000
 C                 -1.01440000    5.98310000    1.71140000
 C                 -0.03230000    6.97140000    1.68340000
 C                  1.41810000    6.98780000    1.68700000
 C                  2.42010000    6.02130000    1.71470000
 C                  3.76060000    6.43980000    1.70620000
 C                  4.80770000    5.47720000    1.77570000
 C                  5.70310000    4.65700000    1.88260000
 C                  6.73020000    3.69300000    2.08250000
 C                  6.41570000    2.42620000    2.60900000
 C                  7.41710000    1.49950000    2.86050000
 C                  8.76190000    1.79920000    2.58790000
 C                  9.83500000    0.83180000    2.90270000
 C                 10.92090000    1.21720000    3.66810000
 C                 11.93160000    0.30130000    4.03220000
 C                 11.87060000   -1.01560000    3.63830000
 C                 10.79720000   -1.45230000    2.82750000
 C                 10.68820000   -2.78310000    2.35250000
 C                  9.65830000   -3.13410000    1.52240000
 N                  8.72820000   -2.22340000    1.12210000
 C                  8.78160000   -0.97970000    1.56100000
 C                  9.78260000   -0.52310000    2.44440000
 C                  9.07240000    3.05860000    2.05390000
 C                  8.07610000    3.99250000    1.80290000
 C                  4.07670000    7.81970000    1.66310000
 C                  3.08670000    8.79290000    1.63810000
 C                  1.77150000    8.35000000    1.64970000
 O                  0.67160000    9.14920000    1.61920000
 C                 -0.41220000    8.32660000    1.64100000
 C                 -1.73580000    8.74160000    1.62140000
 C                 -2.70530000    7.74800000    1.64770000
 C                 -5.39020000    2.55290000    2.64660000
 C                 -6.52770000    1.80420000    2.92130000
 C                -10.01790000    2.23440000    3.68030000
 C                -11.17190000    1.57850000    4.16090000
 C                -11.35490000    0.22980000    3.96140000
 C                -10.39130000   -0.50320000    3.23060000
 H                -11.39380000   -2.43580000    3.33460000
 H                 -9.69590000   -3.58680000    1.92560000
 H                 -7.48790000   -0.10390000    1.45900000
 H                 -8.88230000    3.87820000    1.60190000
 H                 -6.86070000    5.21710000    1.11160000
 H                 -0.75110000    4.92460000    1.74270000
 H                  2.17990000    4.95760000    1.74140000
 H                  5.37710000    2.18510000    2.84220000
 H                  7.15350000    0.53830000    3.30830000
 H                 10.98010000    2.24680000    4.02720000
 H                 12.76250000    0.64730000    4.65070000
 H                 12.64230000   -1.72880000    3.93450000
 H                 11.43250000   -3.52880000    2.63590000
 H                  9.55180000   -4.14660000    1.13040000
 H                  8.01870000   -0.29300000    1.19060000
 H                 10.11170000    3.30410000    1.82400000
 H                  8.33470000    4.96850000    1.38840000
 H                  5.12650000    8.11780000    1.65470000
 H                  3.32640000    9.85620000    1.60220000
 H                 -1.99830000    9.79920000    1.57750000
 H                 -3.76030000    8.02440000    1.62750000
 H                 -4.40970000    2.18600000    2.95570000
 H                 -6.42480000    0.86040000    3.46230000
 H                 -9.88460000    3.29770000    3.89040000
 H                -11.91680000    2.15430000    4.71440000
 H                -12.23730000   -0.28320000    4.34860000

 1 4 1.0 43 1.0 83 1.0 158 1.0
 2 20 1.0 59 1.0 112 1.0 187 1.0
 3 4 1.5 27 2.0 28 1.0
 4 5 1.5
 5 6 2.0 29 1.0
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
 21 22 2.0 37 1.0
 22 38 1.0
 23 24 1.5 39 1.0
 24 25 1.0
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
 44 45 2.0 68 1.0
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
 57 58 2.0 74 1.0
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
 81 82 2.0 129 1.5 130 1.0
 82 83 1.5 131 1.0
 83 84 2.0
 84 85 1.5 132 1.0
 85 86 1.5 129 1.5
 86 87 1.0 126 2.0
 87 88 1.5 125 1.5
 88 89 1.5 133 1.0
 89 90 1.5 134 1.0
 90 91 1.5 124 1.5
 91 92 3.0
 92 93 1.5
 93 94 1.5 123 1.5
 94 95 1.5 135 1.0
 95 96 1.0 121 1.5
 96 97 1.5 119 1.5
 97 98 1.5 136 1.0
 98 99 1.5 117 1.5
 99 100 3.0
 100 101 1.5
 101 102 1.5 116 1.5
 102 103 1.5 137 1.0
 103 104 1.5 138 1.0
 104 105 1.0 115 1.5
 105 106 2.0 114 1.5
 106 107 1.5 139 1.0
 107 108 2.0 140 1.0
 108 109 1.5 141 1.0
 109 110 1.5 114 1.5
 110 111 2.0 142 1.0
 111 112 1.5 143 1.0
 112 113 2.0
 113 114 1.5 144 1.0
 114
 115 116 1.5 145 1.0
 116 146 1.0
 117 118 1.5 147 1.0
 118 119 1.5 148 1.0
 119 120 1.0
 120 121 1.0
 121 122 1.5
 122 123 1.5 149 1.0
 123 150 1.0
 124 125 1.5 151 1.0
 125 152 1.0
 126 127 1.5 153 1.0
 127 128 2.0 154 1.0
 128 129 1.5 155 1.0
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
 146
 147
 148
 149
 150
 151
 152
 153
 154
 155
 156 157 2.0 204 1.5 205 1.0
 157 158 1.5 206 1.0
 158 159 2.0
 159 160 1.5 207 1.0
 160 161 1.5 204 1.5
 161 162 1.0 201 2.0
 162 163 1.5 200 1.5
 163 164 1.5 208 1.0
 164 165 1.5 209 1.0
 165 166 1.5 199 1.5
 166 167 3.0
 167 168 1.5
 168 169 1.5 198 1.5
 169 170 1.5 210 1.0
 170 171 1.0 196 1.5
 171 172 1.5 194 1.5
 172 173 1.5 211 1.0
 173 174 1.5 192 1.5
 174 175 3.0
 175 176 1.5
 176 177 1.5 191 1.5
 177 178 1.5 212 1.0
 178 179 1.5 213 1.0
 179 180 1.0 190 1.5
 180 181 2.0 189 1.5
 181 182 1.5 214 1.0
 182 183 2.0 215 1.0
 183 184 1.5 216 1.0
 184 185 1.5 189 1.5
 185 186 2.0 217 1.0
 186 187 1.5 218 1.0
 187 188 2.0
 188 189 1.5 219 1.0
 189
 190 191 1.5 220 1.0
 191 221 1.0
 192 193 1.5 222 1.0
 193 194 1.5 223 1.0
 194 195 1.0
 195 196 1.0
 196 197 1.5
 197 198 1.5 224 1.0
 198 225 1.0
 199 200 1.5 226 1.0
 200 227 1.0
 201 202 1.5 228 1.0
 202 203 2.0 229 1.0
 203 204 1.5 230 1.0
 204
 205
 206
 207
 208
 209
 210
 211
 212
 213
 214
 215
 216
 217
 218
 219
 220
 221
 222
 223
 224
 225
 226
 227
 228
 229
 230

