%chk=m3_lb_FINAL.chk
%nprocshared=48
%mem=120GB
# SP PBE1PBE/def2svp empiricaldispersion=GD3BJ scrf=(solvent=dmso)

m3_lb

6 1
 Pd                -7.08910000   -4.23760000    0.64250000
 Pd                -5.28410000    4.55010000   -0.22210000
 Pd                11.57010000   -1.41450000   -2.01090000
 C                 -6.96390000   -5.81650000   -3.29260000
 C                 -7.35420000   -5.23310000   -2.11850000
 N                 -6.45550000   -4.96990000   -1.13080000
 C                 -5.17220000   -5.23840000   -1.29430000
 C                 -4.65980000   -5.79430000   -2.48650000
 C                 -3.26500000   -6.05030000   -2.68500000
 C                 -2.22810000   -5.64380000   -1.71160000
 C                 -2.16060000   -4.33400000   -1.20800000
 C                 -1.13670000   -3.95260000   -0.35010000
 C                 -0.15230000   -4.88050000    0.03490000
 C                  0.87280000   -4.52500000    0.95930000
 C                  1.72490000   -4.25700000    1.78670000
 C                  2.73040000   -3.95120000    2.75970000
 N                  3.54760000   -2.93450000    2.47210000
 C                  4.50060000   -2.60700000    3.35220000
 C                  5.35740000   -1.53160000    2.95130000
 C                  6.05610000   -0.65110000    2.48320000
 C                  6.87300000    0.31290000    1.82650000
 C                  8.21940000    0.01620000    1.55030000
 C                  8.97780000    0.87260000    0.76760000
 C                  8.42360000    2.04390000    0.21890000
 C                  9.17260000    2.91040000   -0.71870000
 C                  9.06800000    4.28910000   -0.65020000
 C                  9.81640000    5.13330000   -1.49850000
 C                 10.70910000    4.62120000   -2.41290000
 C                 10.84750000    3.21890000   -2.53790000
 C                 11.76340000    2.60800000   -3.43410000
 C                 11.88730000    1.24300000   -3.47680000
 N                 11.08470000    0.44230000   -2.72530000
 C                 10.17110000    0.98190000   -1.93460000
 C                 10.04430000    2.36610000   -1.72120000
 C                  7.11490000    2.38170000    0.59370000
 C                  6.34180000    1.53280000    1.37390000
 C                  4.67370000   -3.28320000    4.56910000
 C                  3.82250000   -4.34430000    4.85880000
 C                  2.83260000   -4.69340000    3.94740000
 C                 -0.21190000   -6.18970000   -0.47340000
 C                 -1.23220000   -6.55890000   -1.33800000
 C                 -2.87870000   -6.68050000   -3.85560000
 C                 -3.81190000   -7.03100000   -4.85370000
 C                 -5.14960000   -6.74950000   -4.69920000
 C                 -5.59970000   -6.12980000   -3.51010000
 H                 -7.70470000   -6.04040000   -4.06220000
 H                 -8.39230000   -4.96480000   -1.91580000
 H                 -4.51650000   -5.02500000   -0.44780000
 H                 -2.89950000   -3.59250000   -1.51990000
 H                 -1.08950000   -2.92830000    0.02610000
 H                  8.66250000   -0.90350000    1.93830000
 H                 10.02440000    0.62350000    0.59370000
 H                  8.42540000    4.74200000    0.10680000
 H                  9.70340000    6.21530000   -1.39570000
 H                 11.31660000    5.27730000   -3.04000000
 H                 12.39600000    3.22970000   -4.07050000
 H                 12.62770000    0.74270000   -4.10300000
 H                  9.52560000    0.29380000   -1.39160000
 H                  6.67520000    3.30830000    0.22810000
 H                  5.30860000    1.79890000    1.60720000
 H                  5.46100000   -2.97660000    5.25920000
 H                  3.92830000   -4.89590000    5.79560000
 H                  2.14430000   -5.51760000    4.13940000
 H                  0.54780000   -6.91640000   -0.17920000
 H                 -1.26950000   -7.58130000   -1.72060000
 H                 -1.81830000   -6.87660000   -4.02800000
 H                 -3.45520000   -7.51230000   -5.76740000
 H                 -5.87430000   -6.99590000   -5.47800000
 C                 -9.30630000   -0.80230000   -0.69710000
 C                 -8.84580000   -1.89640000   -0.01630000
 N                 -7.65890000   -2.49540000   -0.31550000
 C                 -6.93360000   -2.02320000   -1.31850000
 C                 -7.33330000   -0.92090000   -2.10760000
 C                 -6.53360000   -0.42040000   -3.18380000
 C                 -5.28470000   -1.10160000   -3.59140000
 C                 -4.06600000   -0.41190000   -3.54760000
 C                 -2.87420000   -1.05140000   -3.86360000
 C                 -2.88100000   -2.40120000   -4.25110000
 C                 -1.66330000   -3.12120000   -4.41120000
 C                 -0.64720000   -3.79170000   -4.42490000
 C                  0.50090000   -4.64550000   -4.40800000
 N                  1.44630000   -4.38700000   -3.49480000
 C                  2.51470000   -5.19130000   -3.45220000
 C                  3.54740000   -4.91450000   -2.50100000
 C                  4.51480000   -4.76670000   -1.77510000
 C                  5.71090000   -4.69130000   -1.00580000
 C                  6.88490000   -5.26450000   -1.52980000
 C                  8.05490000   -5.26580000   -0.78620000
 C                  8.08230000   -4.72360000    0.50960000
 C                  9.31520000   -4.84750000    1.31970000
 C                  9.30470000   -5.50790000    2.53380000
 C                 10.50320000   -5.76410000    3.23790000
 C                 11.72490000   -5.38030000    2.73200000
 C                 11.77970000   -4.65960000    1.51600000
 C                 12.98570000   -4.19390000    0.93450000
 C                 12.95270000   -3.41410000   -0.19050000
 N                 11.77750000   -3.03660000   -0.76910000
 C                 10.63880000   -3.51060000   -0.29140000
 C                 10.56720000   -4.36100000    0.82960000
 C                  6.91710000   -4.13690000    1.02110000
 C                  5.74760000   -4.10190000    0.26880000
 C                  2.68030000   -6.28940000   -4.31700000
 C                  1.69160000   -6.54950000   -5.25700000
 C                  0.57930000   -5.71840000   -5.31210000
 C                 -4.11080000   -3.07180000   -4.37450000
 C                 -5.29460000   -2.43410000   -4.03110000
 C                 -6.94880000    0.72650000   -3.83290000
 C                 -8.15490000    1.37620000   -3.48600000
 C                 -8.95810000    0.88780000   -2.48190000
 C                 -8.55410000   -0.26500000   -1.76730000
 H                -10.25660000   -0.34890000   -0.40960000
 H                 -9.42460000   -2.31730000    0.80150000
 H                 -5.97790000   -2.51100000   -1.51150000
 H                 -4.05300000    0.63000000   -3.22320000
 H                 -1.92560000   -0.51760000   -3.78670000
 H                  6.86470000   -5.72150000   -2.52050000
 H                  8.95390000   -5.72910000   -1.19980000
 H                  8.35790000   -5.87940000    2.93290000
 H                 10.45310000   -6.30480000    4.18610000
 H                 12.65280000   -5.61240000    3.25940000
 H                 13.94590000   -4.44030000    1.39190000
 H                 13.86940000   -3.04530000   -0.64940000
 H                  9.71790000   -3.17970000   -0.77130000
 H                  6.93040000   -3.69680000    2.02080000
 H                  4.85010000   -3.62970000    0.67480000
 H                  3.57130000   -6.91410000   -4.23940000
 H                  1.78770000   -7.39360000   -5.94360000
 H                 -0.22070000   -5.88120000   -6.03600000
 H                 -4.12440000   -4.10480000   -4.72290000
 H                 -6.24100000   -2.97620000   -4.10460000
 H                 -6.34550000    1.11920000   -4.65380000
 H                 -8.45680000    2.26770000   -4.03990000
 H                 -9.90770000    1.36500000   -2.22990000
 C                 -7.37780000   -8.30570000    1.71370000
 C                 -7.53060000   -7.06450000    1.16010000
 N                 -6.66800000   -6.04980000    1.44740000
 C                 -5.64250000   -6.25070000    2.25210000
 C                 -5.39600000   -7.49290000    2.87900000
 C                 -4.32060000   -7.69010000    3.80490000
 C                 -3.34810000   -6.61950000    4.11850000
 C                 -2.64410000   -5.95210000    3.10350000
 C                 -1.72750000   -4.96030000    3.41620000
 C                 -1.47880000   -4.60800000    4.75700000
 C                 -0.53260000   -3.57300000    5.00510000
 C                  0.29270000   -2.67760000    5.05360000
 C                  1.23450000   -1.62000000    4.84080000
 N                  1.21990000   -1.13740000    3.59370000
 C                  2.08210000   -0.17340000    3.25690000
 C                  1.92350000    0.29870000    1.91690000
 C                  1.53740000    0.65980000    0.81930000
 C                  0.83880000    1.08560000   -0.34450000
 C                 -0.45640000    0.57220000   -0.54600000
 C                 -1.25680000    1.07110000   -1.56040000
 C                 -0.76950000    2.04910000   -2.44280000
 C                 -1.59340000    2.46740000   -3.59970000
 C                 -1.12250000    2.24770000   -4.88120000
 C                 -1.92620000    2.49190000   -6.01590000
 C                 -3.21730000    2.94630000   -5.88120000
 C                 -3.73040000    3.22000000   -4.59200000
 C                 -5.04940000    3.68980000   -4.38140000
 C                 -5.47850000    3.97600000   -3.11480000
 N                 -4.65530000    3.86140000   -2.03630000
 C                 -3.43070000    3.38590000   -2.18470000
 C                 -2.90790000    3.01020000   -3.44330000
 C                  0.53190000    2.53610000   -2.25890000
 C                  1.32220000    2.07900000   -1.20970000
 C                  3.00570000    0.36370000    4.16640000
 C                  3.00630000   -0.13210000    5.46620000
 C                  2.11480000   -1.14050000    5.82280000
 C                 -2.16970000   -5.28320000    5.77620000
 C                 -3.09010000   -6.27390000    5.45450000
 C                 -4.21780000   -8.91820000    4.43360000
 C                 -5.11630000   -9.97060000    4.15630000
 C                 -6.14240000   -9.80460000    3.25570000
 C                 -6.31010000   -8.55700000    2.61080000
 H                 -8.08950000   -9.09830000    1.47570000
 H                 -8.34710000   -6.82940000    0.47550000
 H                 -4.98540000   -5.40120000    2.43520000
 H                 -2.79670000   -6.22760000    2.05680000
 H                 -1.18130000   -4.45010000    2.62250000
 H                 -0.84050000   -0.19530000    0.12840000
 H                 -2.26920000    0.68330000   -1.68290000
 H                 -0.12110000    1.83180000   -5.01500000
 H                 -1.51860000    2.28730000   -7.00880000
 H                 -3.85610000    3.10520000   -6.75270000
 H                 -5.72710000    3.82050000   -5.22710000
 H                 -6.49490000    4.31920000   -2.92360000
 H                 -2.81160000    3.31530000   -1.28740000
 H                  0.92170000    3.30190000   -2.93390000
 H                  2.31120000    2.50830000   -1.04480000
 H                  3.68790000    1.15650000    3.85700000
 H                  3.70730000    0.26550000    6.20350000
 H                  2.09630000   -1.55650000    6.83120000
 H                 -1.98820000   -5.02190000    6.82060000
 H                 -3.63400000   -6.78230000    6.25410000
 H                 -3.40640000   -9.08490000    5.14560000
 H                 -4.98420000  -10.92850000    4.66490000
 H                 -6.83990000  -10.61570000    3.03590000
 C                 -7.21210000   -3.43270000    4.83560000
 C                 -6.79810000   -3.52680000    3.53270000
 N                 -7.68790000   -3.66850000    2.50940000
 C                 -8.98080000   -3.63720000    2.76210000
 C                 -9.50970000   -3.43390000    4.05290000
 C                -10.90870000   -3.22150000    4.26330000
 C                -11.72810000   -2.91030000    3.06170000
 C                -11.91210000   -1.55220000    2.76310000
 C                -12.38890000   -1.15730000    1.52280000
 C                -12.69650000   -2.12000000    0.54370000
 C                -12.97660000   -1.66610000   -0.77750000
 C                -13.08290000   -1.16550000   -1.88310000
 C                -13.10430000   -0.48800000   -3.14460000
 N                -12.51430000    0.71440000   -3.15990000
 C                -12.47840000    1.39660000   -4.30980000
 C                -11.81960000    2.67060000   -4.25880000
 C                -11.23110000    3.73010000   -4.13900000
 C                -10.49570000    4.93390000   -3.92120000
 C                -11.01570000    5.95310000   -3.10420000
 C                -10.24730000    7.07210000   -2.80140000
 C                 -8.94090000    7.19170000   -3.29010000
 C                 -8.05750000    8.30750000   -2.87550000
 C                 -8.31470000    9.62150000   -3.20990000
 C                 -7.40980000   10.65180000   -2.86370000
 C                 -6.23930000   10.37970000   -2.19230000
 C                 -5.94900000    9.05060000   -1.80200000
 C                 -4.77700000    8.68380000   -1.09240000
 C                 -4.58530000    7.38280000   -0.70880000
 N                 -5.51650000    6.42220000   -0.96470000
 C                 -6.61260000    6.72070000   -1.63350000
 C                 -6.87710000    8.01710000   -2.12320000
 C                 -8.44180000    6.19640000   -4.14440000
 C                 -9.20860000    5.08570000   -4.46690000
 C                -13.03580000    0.90580000   -5.50110000
 C                -13.64940000   -0.34210000   -5.47590000
 C                -13.69010000   -1.05930000   -4.28550000
 C                -12.60570000   -3.48300000    0.87440000
 C                -12.12400000   -3.87040000    2.12210000
 C                -11.36210000   -3.09430000    5.55840000
 C                -10.46530000   -3.16220000    6.65250000
 C                 -9.10790000   -3.28350000    6.45980000
 C                 -8.59630000   -3.39810000    5.14440000
 H                 -6.47010000   -3.36110000    5.63340000
 H                 -5.74130000   -3.50300000    3.26600000
 H                 -9.65190000   -3.80290000    1.91940000
 H                -11.61650000   -0.79680000    3.49510000
 H                -12.48180000   -0.09610000    1.28280000
 H                -12.01960000    5.84950000   -2.68740000
 H                -10.65110000    7.84480000   -2.14280000
 H                 -9.21590000    9.86540000   -3.77740000
 H                 -7.64170000   11.67880000   -3.15630000
 H                 -5.53120000   11.17480000   -1.94870000
 H                 -4.02610000    9.43730000   -0.84690000
 H                 -3.69440000    7.06360000   -0.16480000
 H                 -7.33930000    5.91950000   -1.77810000
 H                 -7.43170000    6.29120000   -4.55040000
 H                 -8.80580000    4.31540000   -5.12770000
 H                -12.98350000    1.49770000   -6.41630000
 H                -14.09510000   -0.75540000   -6.38350000
 H                -14.16120000   -2.04170000   -4.22680000
 H                -12.86770000   -4.23850000    0.13090000
 H                -12.00150000   -4.93270000    2.34740000
 H                -12.42230000   -2.90700000    5.74360000
 H                -10.86290000   -3.07280000    7.66640000
 H                 -8.41560000   -3.27590000    7.30480000
 C                 -5.74880000    0.45400000    0.72550000
 C                 -5.94140000    1.73360000    0.28550000
 N                 -5.03800000    2.71060000    0.56390000
 C                 -3.95850000    2.45150000    1.28010000
 C                 -3.62720000    1.14810000    1.71290000
 C                 -2.41220000    0.83990000    2.40730000
 C                 -1.35890000    1.83300000    2.72280000
 C                 -0.87800000    2.75430000    1.77750000
 C                  0.26300000    3.50420000    2.02850000
 C                  0.95290000    3.36340000    3.24420000
 C                  2.20050000    4.02730000    3.43430000
 C                  3.27890000    4.58370000    3.53980000
 C                  4.48680000    5.34960000    3.64520000
 N                  5.64680000    4.74330000    3.36330000
 C                  6.76980000    5.47070000    3.43230000
 C                  8.01570000    4.82710000    3.13780000
 C                  9.11340000    4.36170000    2.88580000
 C                 10.38870000    3.85860000    2.49690000
 C                 10.77430000    2.53620000    2.77260000
 C                 11.98050000    2.04110000    2.29090000
 C                 12.83480000    2.85110000    1.52890000
 C                 14.12950000    2.36170000    0.99960000
 C                 15.28900000    3.04550000    1.31920000
 C                 16.54770000    2.65040000    0.81840000
 C                 16.66040000    1.57720000   -0.03340000
 C                 15.50060000    0.86280000   -0.41340000
 C                 15.54700000   -0.21020000   -1.33660000
 C                 14.39340000   -0.84250000   -1.70760000
 N                 13.18400000   -0.47960000   -1.19400000
 C                 13.09740000    0.51230000   -0.32650000
 C                 14.22680000    1.24060000    0.11340000
 C                 12.46460000    4.18370000    1.28930000
 C                 11.26130000    4.68370000    1.76400000
 C                  6.77630000    6.83490000    3.77270000
 C                  5.56830000    7.45090000    4.07160000
 C                  4.39940000    6.70290000    4.01350000
 C                  0.42180000    2.51910000    4.23320000
 C                 -0.71270000    1.76560000    3.96910000
 C                 -2.18500000   -0.48410000    2.75480000
 C                 -3.11410000   -1.50530000    2.46330000
 C                 -4.29100000   -1.22350000    1.80980000
 C                 -4.56390000    0.10800000    1.41840000
 H                 -6.50630000   -0.30420000    0.53920000
 H                 -6.81870000    2.02690000   -0.29440000
 H                 -3.32040000    3.30290000    1.52100000
 H                 -1.34240000    2.82980000    0.79380000
 H                  0.66580000    4.16240000    1.25640000
 H                 10.11300000    1.89300000    3.35540000
 H                 12.26530000    1.01020000    2.51480000
 H                 15.22870000    3.90680000    1.98840000
 H                 17.43720000    3.21230000    1.11290000
 H                 17.62880000    1.27050000   -0.43430000
 H                 16.50060000   -0.52080000   -1.76720000
 H                 14.38740000   -1.65960000   -2.43180000
 H                 12.09990000    0.77330000    0.02880000
 H                 13.12150000    4.82750000    0.70000000
 H                 10.97400000    5.71450000    1.54850000
 H                  7.71890000    7.38350000    3.80420000
 H                  5.53720000    8.50810000    4.34460000
 H                  3.42690000    7.14510000    4.23510000
 H                  0.92800000    2.42660000    5.19560000
 H                 -1.09380000    1.08730000    4.73550000
 H                 -1.23380000   -0.75590000    3.22000000
 H                 -2.87820000   -2.53210000    2.75230000
 H                 -5.01360000   -2.00780000    1.57040000
 C                 -7.49280000    5.65110000    3.25790000
 C                 -7.13790000    5.17770000    2.02230000
 N                 -5.85270000    5.25340000    1.58240000
 C                 -4.92030000    5.81330000    2.32760000
 C                 -5.19020000    6.36980000    3.59520000
 C                 -4.16800000    6.98900000    4.38200000
 C                 -2.80690000    7.14180000    3.82040000
 C                 -1.70070000    6.55730000    4.45060000
 C                 -0.44730000    6.58760000    3.85120000
 C                 -0.27360000    7.20860000    2.60180000
 C                  0.97680000    7.11390000    1.92250000
 C                  2.02410000    6.93960000    1.32730000
 C                  3.25890000    6.66880000    0.65180000
 N                  3.35820000    5.48300000    0.03770000
 C                  4.50700000    5.18740000   -0.58210000
 C                  4.63050000    3.90780000   -1.20960000
 C                  4.86720000    2.81410000   -1.69180000
 C                  5.24530000    1.55040000   -2.22650000
 C                  6.33480000    1.48070000   -3.11430000
 C                  6.74860000    0.25750000   -3.62050000
 C                  6.08990000   -0.93010000   -3.26140000
 C                  6.49960000   -2.22320000   -3.85410000
 C                  5.56800000   -3.01790000   -4.49590000
 C                  5.94010000   -4.21820000   -5.14040000
 C                  7.25020000   -4.63760000   -5.15660000
 C                  8.23570000   -3.87460000   -4.48670000
 C                  9.59670000   -4.26360000   -4.42640000
 C                 10.49680000   -3.51170000   -3.72110000
 N                 10.11420000   -2.39030000   -3.04580000
 C                  8.85670000   -1.99150000   -3.08720000
 C                  7.85880000   -2.67270000   -3.81620000
 C                  5.00480000   -0.85770000   -2.37450000
 C                  4.58820000    0.36290000   -1.85640000
 C                  5.60500000    6.06620000   -0.61160000
 C                  5.49310000    7.29030000    0.03500000
 C                  4.30250000    7.60710000    0.67850000
 C                 -1.36670000    7.85580000    2.00070000
 C                 -2.61460000    7.82120000    2.60680000
 C                 -4.49020000    7.41260000    5.65630000
 C                 -5.80510000    7.28380000    6.15880000
 C                 -6.81000000    6.73170000    5.39710000
 C                 -6.52090000    6.25500000    4.09620000
 H                 -8.52600000    5.55920000    3.59810000
 H                 -7.85690000    4.71200000    1.34620000
 H                 -3.90820000    5.81770000    1.92020000
 H                 -1.83380000    6.03370000    5.40030000
 H                  0.40190000    6.09960000    4.33290000
 H                  6.84960000    2.39770000   -3.40810000
 H                  7.58110000    0.22620000   -4.32750000
 H                  4.52580000   -2.69390000   -4.53290000
 H                  5.17010000   -4.80810000   -5.64340000
 H                  7.54420000   -5.55590000   -5.66950000
 H                  9.92980000   -5.16910000   -4.93690000
 H                 11.55110000   -3.78380000   -3.65120000
 H                  8.59200000   -1.10800000   -2.50720000
 H                  4.48880000   -1.77420000   -2.07720000
 H                  3.75480000    0.40050000   -1.15230000
 H                  6.51940000    5.77750000   -1.13170000
 H                  6.32960000    7.99280000    0.03850000
 H                  4.16990000    8.55590000    1.20020000
 H                 -1.23500000    8.36490000    1.04370000
 H                 -3.46170000    8.31230000    2.12170000
 H                 -3.72140000    7.87950000    6.27640000
 H                 -6.02320000    7.64460000    7.16690000
 H                 -7.82880000    6.64850000    5.78180000

 1 6 1.0 71 1.0 136 1.0 201 1.0
 2 162 1.0 227 1.0 266 1.0 331 1.0
 3 32 1.0 97 1.0 292 1.0 357 1.0
 4 5 2.0 45 1.5 46 1.0
 5 6 1.5 47 1.0
 6 7 2.0
 7 8 1.5 48 1.0
 8 9 1.5 45 1.5
 9 10 1.0 42 2.0
 10 11 1.5 41 1.5
 11 12 1.5 49 1.0
 12 13 1.5 50 1.0
 13 14 1.5 40 1.5
 14 15 3.0
 15 16 1.5
 16 17 1.5 39 1.5
 17 18 1.5
 18 19 1.5 37 1.5
 19 20 3.0
 20 21 1.5
 21 22 1.5 36 1.5
 22 23 1.5 51 1.0
 23 24 1.5 52 1.0
 24 25 1.0 35 1.5
 25 26 2.0 34 1.5
 26 27 1.5 53 1.0
 27 28 2.0 54 1.0
 28 29 1.5 55 1.0
 29 30 1.5 34 1.5
 30 31 2.0 56 1.0
 31 32 1.5 57 1.0
 32 33 1.5
 33 34 1.5 58 1.0
 34
 35 36 1.5 59 1.0
 36 60 1.0
 37 38 1.5 61 1.0
 38 39 1.5 62 1.0
 39 63 1.0
 40 41 1.5 64 1.0
 41 65 1.0
 42 43 1.5 66 1.0
 43 44 2.0 67 1.0
 44 45 1.5 68 1.0
 45
 46
 47
 48
 49
 50
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
 63
 64
 65
 66
 67
 68
 69 70 2.0 110 1.5 111 1.0
 70 71 1.5 112 1.0
 71 72 1.5
 72 73 1.5 113 1.0
 73 74 1.5 110 1.5
 74 75 1.0 107 2.0
 75 76 1.5 106 1.5
 76 77 1.5 114 1.0
 77 78 1.5 115 1.0
 78 79 1.5 105 1.5
 79 80 3.0
 80 81 1.5
 81 82 1.5 104 1.5
 82 83 1.5
 83 84 1.5 102 1.5
 84 85 3.0
 85 86 1.5
 86 87 1.5 101 1.5
 87 88 1.5 116 1.0
 88 89 1.5 117 1.0
 89 90 1.0 100 1.5
 90 91 2.0 99 1.5
 91 92 1.5 118 1.0
 92 93 2.0 119 1.0
 93 94 1.5 120 1.0
 94 95 1.5 99 1.5
 95 96 2.0 121 1.0
 96 97 1.5 122 1.0
 97 98 2.0
 98 99 1.5 123 1.0
 99
 100 101 1.5 124 1.0
 101 125 1.0
 102 103 1.5 126 1.0
 103 104 1.5 127 1.0
 104 128 1.0
 105 106 1.5 129 1.0
 106 130 1.0
 107 108 1.5 131 1.0
 108 109 2.0 132 1.0
 109 110 1.5 133 1.0
 110
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
 134 135 2.0 175 1.5 176 1.0
 135 136 1.5 177 1.0
 136 137 2.0
 137 138 1.5 178 1.0
 138 139 1.5 175 1.5
 139 140 1.0 172 2.0
 140 141 1.5 171 1.5
 141 142 1.5 179 1.0
 142 143 1.5 180 1.0
 143 144 1.5 170 1.5
 144 145 3.0
 145 146 1.5
 146 147 1.5 169 1.5
 147 148 1.5
 148 149 1.5 167 1.5
 149 150 3.0
 150 151 1.5
 151 152 1.5 166 1.5
 152 153 2.0 181 1.0
 153 154 1.5 182 1.0
 154 155 1.0 165 1.5
 155 156 2.0 164 1.5
 156 157 1.5 183 1.0
 157 158 2.0 184 1.0
 158 159 1.5 185 1.0
 159 160 1.5 164 1.5
 160 161 2.0 186 1.0
 161 162 1.5 187 1.0
 162 163 2.0
 163 164 1.5 188 1.0
 164
 165 166 1.5 189 1.0
 166 190 1.0
 167 168 1.5 191 1.0
 168 169 1.5 192 1.0
 169 193 1.0
 170 171 1.5 194 1.0
 171 195 1.0
 172 173 1.5 196 1.0
 173 174 2.0 197 1.0
 174 175 1.5 198 1.0
 175
 176
 177
 178
 179
 180
 181
 182
 183
 184
 185
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
 199 200 2.0 240 1.5 241 1.0
 200 201 1.5 242 1.0
 201 202 2.0
 202 203 1.5 243 1.0
 203 204 1.5 240 1.5
 204 205 1.0 237 2.0
 205 206 1.5 236 1.5
 206 207 1.5 244 1.0
 207 208 1.5 245 1.0
 208 209 1.5 235 1.5
 209 210 3.0
 210 211 1.5
 211 212 1.5 234 1.5
 212 213 1.5
 213 214 1.5 232 1.5
 214 215 3.0
 215 216 1.5
 216 217 1.5 231 1.5
 217 218 1.5 246 1.0
 218 219 1.5 247 1.0
 219 220 1.0 230 1.5
 220 221 2.0 229 1.5
 221 222 1.5 248 1.0
 222 223 2.0 249 1.0
 223 224 1.5 250 1.0
 224 225 1.5 229 1.5
 225 226 2.0 251 1.0
 226 227 1.5 252 1.0
 227 228 2.0
 228 229 1.5 253 1.0
 229
 230 231 1.5 254 1.0
 231 255 1.0
 232 233 1.5 256 1.0
 233 234 1.5 257 1.0
 234 258 1.0
 235 236 1.5 259 1.0
 236 260 1.0
 237 238 1.5 261 1.0
 238 239 2.0 262 1.0
 239 240 1.5 263 1.0
 240
 241
 242
 243
 244
 245
 246
 247
 248
 249
 250
 251
 252
 253
 254
 255
 256
 257
 258
 259
 260
 261
 262
 263
 264 265 2.0 305 1.5 306 1.0
 265 266 1.5 307 1.0
 266 267 2.0
 267 268 1.5 308 1.0
 268 269 1.5 305 1.5
 269 270 1.0 302 1.5
 270 271 1.5 301 1.5
 271 272 1.5 309 1.0
 272 273 1.5 310 1.0
 273 274 1.5 300 1.5
 274 275 3.0
 275 276 1.5
 276 277 1.5 299 1.5
 277 278 1.5
 278 279 1.5 297 1.5
 279 280 3.0
 280 281 1.5
 281 282 1.5 296 1.5
 282 283 1.5 311 1.0
 283 284 1.5 312 1.0
 284 285 1.0 295 1.5
 285 286 2.0 294 1.5
 286 287 1.5 313 1.0
 287 288 2.0 314 1.0
 288 289 1.5 315 1.0
 289 290 1.5 294 1.5
 290 291 2.0 316 1.0
 291 292 1.5 317 1.0
 292 293 2.0
 293 294 1.5 318 1.0
 294
 295 296 1.5 319 1.0
 296 320 1.0
 297 298 1.5 321 1.0
 298 299 1.5 322 1.0
 299 323 1.0
 300 301 1.5 324 1.0
 301 325 1.0
 302 303 1.5 326 1.0
 303 304 2.0 327 1.0
 304 305 1.5 328 1.0
 305
 306
 307
 308
 309
 310
 311
 312
 313
 314
 315
 316
 317
 318
 319
 320
 321
 322
 323
 324
 325
 326
 327
 328
 329 330 2.0 370 1.5 371 1.0
 330 331 1.5 372 1.0
 331 332 2.0
 332 333 1.5 373 1.0
 333 334 1.5 370 1.5
 334 335 1.0 367 2.0
 335 336 1.5 366 1.5
 336 337 1.5 374 1.0
 337 338 1.5 375 1.0
 338 339 1.5 365 1.5
 339 340 3.0
 340 341 1.5
 341 342 1.5 364 1.5
 342 343 1.5
 343 344 1.5 362 1.5
 344 345 3.0
 345 346 1.5
 346 347 1.5 361 1.5
 347 348 1.5 376 1.0
 348 349 1.5 377 1.0
 349 350 1.0 360 1.5
 350 351 2.0 359 1.5
 351 352 1.5 378 1.0
 352 353 2.0 379 1.0
 353 354 1.5 380 1.0
 354 355 1.5 359 1.5
 355 356 2.0 381 1.0
 356 357 1.5 382 1.0
 357 358 2.0
 358 359 1.5 383 1.0
 359
 360 361 1.5 384 1.0
 361 385 1.0
 362 363 1.5 386 1.0
 363 364 1.5 387 1.0
 364 388 1.0
 365 366 1.5 389 1.0
 366 390 1.0
 367 368 1.5 391 1.0
 368 369 2.0 392 1.0
 369 370 1.5 393 1.0
 370
 371
 372
 373
 374
 375
 376
 377
 378
 379
 380
 381
 382
 383
 384
 385
 386
 387
 388
 389
 390
 391
 392
 393

