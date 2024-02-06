%nprocshared=48
%mem=120GB
%chk=m2_lb_xray.chk
# SP HSEH1PBE/def2svp scrf=(solvent=dmso)

m2_lb

4 1
 C                 -8.10798200    1.73198700   -2.32116900
 C                 -7.93147600    2.33186500   -3.54029900
 C                 -6.88073900    1.91683200   -4.39630100
 C                 -6.59423100    2.54221700   -5.63226900
 C                 -5.51982200    2.11083100   -6.37802600
 C                 -4.72026900    1.02506100   -5.95293800
 C                 -4.98856000    0.35569700   -4.77530200
 C                 -6.06346800    0.83072700   -3.96079400
 C                 -6.31231500    0.29668300   -2.68244900
 C                 -4.19828100   -0.82608600   -4.36007800
 C                 -2.82195400   -0.71820500   -4.12667100
 C                 -2.08425800   -1.82028500   -3.71447200
 C                 -2.71185900   -3.06404000   -3.52574100
 C                 -4.08800200   -3.18191100   -3.78977100
 C                 -4.82009800   -2.07432700   -4.19605800
 C                 -1.97906200   -4.16672400   -2.99959800
 C                 -1.36415000   -5.06208700   -2.44889200
 C                 -0.64688700   -5.98018100   -1.61443300
 C                 -0.66296500   -7.36844700   -1.82151700
 C                  0.02976300   -8.17421100   -0.92252900
 C                  0.71204000   -7.58871900    0.13999800
 C                  0.67535400   -6.19062900    0.26027900
 C                  1.38598500   -5.48088100    1.28229100
 C                  1.99878100   -4.73097000    2.02062000
 C                  2.73024700   -3.77280400    2.77984800
 C                  4.10421100   -3.95083800    3.02084800
 C                  4.83655500   -2.95999300    3.66103500
 C                  4.21681500   -1.77392700    4.08628900
 C                  2.84206800   -1.61325300    3.87515100
 C                  2.10415500   -2.59740500    3.23038900
 C                  5.00843500   -0.70995300    4.74608900
 C                  4.74759900   -0.31495900    6.04324200
 C                  5.54906300    0.65381600    6.68996900
 C                  6.61787500    1.24030300    6.04949000
 C                  6.89648300    0.90066000    4.70498900
 C                  7.94152100    1.49423100    3.95377500
 C                  8.11024800    1.17587600    2.63191700
 C                  6.31761800   -0.30720100    2.68201900
 C                  6.07700100   -0.06567300    4.04785700
 C                 -8.12309900   -2.31324100   -1.72534300
 C                 -7.95427600   -3.53343800   -2.32538100
 C                 -6.90505600   -4.39354900   -1.91468900
 C                 -6.62689400   -5.63121200   -2.54047900
 C                 -5.55379700   -6.38141200   -2.11357400
 C                 -4.74773100   -5.95932400   -1.03146800
 C                 -5.00823700   -4.78058700   -0.36091600
 C                 -6.08126100   -3.96118100   -0.83217800
 C                 -6.32218700   -2.68124600   -0.29811200
 C                 -4.21153500   -4.36913300    0.81781200
 C                 -2.83333700   -4.14826600    0.70617300
 C                 -2.08958600   -3.73721400    1.80465500
 C                 -2.71297400   -3.53772600    3.04884300
 C                 -4.09070800   -3.79124900    3.17120700
 C                 -4.82878800   -4.19603500    2.06714700
 C                 -1.97548100   -3.01130500    4.14828600
 C                 -1.35975000   -2.46029200    5.04288900
 C                 -0.64156400   -1.62651900    5.96087000
 C                 -0.65575300   -1.83566200    7.34888200
 C                  0.03684400   -0.93699300    8.15520800
 C                  0.71741300    0.12686000    7.57028700
 C                  0.67861900    0.24935400    6.17246400
 C                  1.38619800    1.27397400    5.46342200
 C                  1.99588400    2.01524000    4.71384700
 C                  2.72708400    2.77798300    3.75830000
 C                  2.10035700    3.23839800    2.58714300
 C                  2.83895000    3.88587500    1.60523900
 C                  4.21504000    4.08958700    1.76376800
 C                  4.83540400    3.65466400    2.94604300
 C                  4.10241300    3.01226000    3.93482800
 C                  5.00831000    4.75089600    0.70209700
 C                  4.74819100    6.04819200    0.30725900
 C                  5.55183000    6.69555300   -0.65939300
 C                  6.62279400    6.05600500   -1.24310600
 C                  6.90127900    4.71148300   -0.90316800
 C                  7.94929500    3.96111400   -1.49288100
 C                  8.11702900    2.63902400   -1.17495900
 C                  6.31756200    2.68725400    0.29969000
 C                  6.07865600    4.05352100    0.05981900
 C                 -8.12724700   -1.72140000    2.31416400
 C                 -7.95512400   -2.32516900    3.53202900
 C                 -6.90263500   -1.91825700    4.39001300
 C                 -6.62137800   -2.54832500    5.62475900
 C                 -5.54571800   -2.12453700    6.37313100
 C                 -4.74034000   -1.04133100    5.95245400
 C                 -5.00387900   -0.36673900    4.77666200
 C                 -6.07933100   -0.83501200    3.95874200
 C                 -6.32370100   -0.29707500    2.68107000
 C                 -4.20890500    0.81380900    4.36711600
 C                 -4.82812100    2.06286500    4.19876300
 C                 -4.09231800    3.16914900    3.79598600
 C                 -2.71483500    3.04933600    3.53965400
 C                 -2.08945200    1.80534800    3.73413600
 C                 -2.83094600    0.70457400    4.14326900
 C                 -1.98008000    4.15123000    3.01459800
 C                 -1.36755400    5.04790000    2.46330900
 C                 -0.65589300    5.96832800    1.62667700
 C                 -0.68171900    7.35704100    1.83023100
 C                  0.00390900    8.16544100    0.92819900
 C                  0.68891900    7.58200000   -0.13366800
 C                  0.66136500    6.18348500   -0.25089100
 C                  1.37273100    5.47593900   -1.27406000
 C                  1.98481400    4.72708500   -2.01405500
 C                  2.71660500    3.77060900   -2.77518600
 C                  2.08815400    2.60154000   -3.23867100
 C                  2.82589800    1.61846100   -3.88521000
 C                  4.20306000    1.77334600   -4.08475000
 C                  4.82542100    2.95320000   -3.64595600
 C                  4.09316000    3.94345900   -3.00490300
 C                  4.99412900    0.71098500   -4.74740800
 C                  4.72966500    0.31623200   -6.04397200
 C                  5.53053500   -0.65072500   -6.69415000
 C                  6.60344100   -1.23488700   -6.05826000
 C                  6.88656600   -0.89499000   -4.71475300
 C                  7.93738100   -1.48499300   -3.96859200
 C                  8.11100100   -1.16683500   -2.64730900
 C                  6.31068100    0.30729400   -2.68811000
 C                  6.06653900    0.06815600   -4.05366700
 C                 -8.13130000    2.31670900    1.72360600
 C                 -7.96042600    3.53273200    2.33142500
 C                 -6.90767200    4.39225500    1.92847800
 C                 -6.62648300    5.62483400    2.56289800
 C                 -5.55128100    6.37515300    2.14140300
 C                 -4.74638300    5.95903100    1.05601900
 C                 -5.00956200    4.78535200    0.37787200
 C                 -6.08396500    3.96503900    0.84409200
 C                 -6.32588100    2.68836200    0.30287600
 C                 -4.21562600    4.37859200   -0.80439700
 C                 -4.83518400    4.21972900   -2.05451200
 C                 -4.10156300    3.81493700   -3.16155900
 C                 -2.72634100    3.54686400   -3.04127500
 C                 -2.10024400    3.73362100   -1.79645500
 C                 -2.83946600    4.14471600   -0.69491200
 C                 -1.99458300    3.01760400   -4.14306500
 C                 -1.38458900    2.46336600   -5.03962100
 C                 -0.67681800    1.62446700   -5.96077900
 C                 -0.71510900    1.81977200   -7.35034500
 C                 -0.03370000    0.91525200   -8.15944900
 C                  0.65911000   -0.14151700   -7.57588500
 C                  0.64343600   -0.25103600   -6.17655400
 C                  1.36034500   -1.27034800   -5.46900700
 C                  1.97483900   -2.00823700   -4.72003400
 C                  2.70903300   -2.76803300   -3.76430100
 C                  2.08403600   -3.22901700   -2.59240600
 C                  2.82416400   -3.87540800   -1.61095100
 C                  4.20036300   -4.07767600   -1.77071800
 C                  4.81926000   -3.64110000   -2.95318500
 C                  4.08463700   -2.99977800   -3.94147400
 C                  4.99463600   -4.74135900   -0.71129500
 C                  4.73242800   -6.03868300   -0.31770100
 C                  5.53738400   -6.68973700    0.64521000
 C                  6.61233400   -6.05424000    1.22602300
 C                  6.89310100   -4.70984500    0.88755500
 C                  7.94545900   -3.96352700    1.47466700
 C                  8.11613100   -2.64137000    1.15832900
 C                  6.31107400   -2.68197600   -0.30987900
 C                  6.06888600   -4.04790900   -0.07133900
 H                 -8.90059900    2.03692000   -1.63746200
 H                 -8.59809600    3.14156300   -3.84068400
 H                 -7.21912500    3.36977800   -5.97299000
 H                 -5.28564400    2.60319200   -7.32415300
 H                 -3.89467800    0.68724600   -6.58319700
 H                 -5.67364400   -0.49918900   -2.29581800
 H                 -2.32311300    0.24422800   -4.24216300
 H                 -1.01395900   -1.71061000   -3.52980900
 H                 -4.57804100   -4.14850500   -3.66432200
 H                 -5.88846800   -2.17949200   -4.40179100
 H                 -1.20562500   -7.79391100   -2.66659200
 H                  0.03834000   -9.25835600   -1.05027800
 H                  1.26299900   -8.19038100    0.86397500
 H                  4.59247500   -4.86933100    2.69190800
 H                  5.90344700   -3.11120800    3.84507500
 H                  2.34546400   -0.69643600    4.19352200
 H                  1.03523300   -2.44676600    3.06782400
 H                  3.92691500   -0.78442000    6.59016100
 H                  5.32112000    0.92667300    7.72243000
 H                  7.24403600    1.97458200    6.55947700
 H                  8.61042500    2.21954200    4.41933700
 H                  8.89853200    1.62441700    2.02674500
 H                  5.67718500   -1.00121400    2.13516200
 H                 -8.91481800   -1.62669500   -2.02655900
 H                 -8.62604800   -3.83205500   -3.13149000
 H                 -7.25761000   -5.97018800   -3.36427800
 H                 -5.32564000   -7.32894900   -2.60606600
 H                 -3.92358000   -6.59315400   -0.69704600
 H                 -5.67846300   -2.29663800    0.49468800
 H                 -2.33789100   -4.27204400   -0.25700500
 H                 -1.01803700   -3.56186800    1.69146300
 H                 -4.57719000   -3.65863200    4.13860800
 H                 -5.89834800   -4.39334600    2.17641600
 H                 -1.19699700   -2.68199500    7.77366300
 H                  0.04644000   -1.06556900    9.23926100
 H                  1.26836300    0.85044100    8.17242400
 H                  1.03041600    3.08149600    2.43767000
 H                  2.34151700    4.21178100    0.69154300
 H                  5.90355700    3.83251400    3.09553600
 H                  4.59100300    2.67571800    4.85032300
 H                  3.92548600    6.59415000    0.77418400
 H                  5.32381200    7.72790300   -0.93265900
 H                  7.25076400    6.56668100   -1.97534700
 H                  8.62146600    4.42760400   -2.21457500
 H                  8.90837500    2.03497100   -1.61979600
 H                  5.67359200    2.13920100    0.98945600
 H                 -8.92164900   -2.01964200    1.62933300
 H                 -8.62625400   -3.13176000    3.83053000
 H                 -7.25168800   -3.37293700    5.96257800
 H                 -5.31493800   -2.62045400    7.31823200
 H                 -3.91439000   -0.70934400    6.58526000
 H                 -5.68021500    0.49606200    2.29673900
 H                 -5.89754700    2.17003300    4.39783000
 H                 -4.58041000    4.13621900    3.66685500
 H                 -1.01812200    1.69415300    3.55621100
 H                 -2.33387000   -0.25827600    4.26312700
 H                 -1.22636200    7.78073800    2.67492800
 H                  0.00466800    9.25002600    1.05260000
 H                  1.23465100    8.18579900   -0.85979300
 H                  1.01740100    2.45463800   -3.08486600
 H                  2.32680300    0.70659300   -4.21366800
 H                  5.89451600    3.10038500   -3.82034600
 H                  4.58342900    4.85730100   -2.66629300
 H                  3.90624700    0.78471700   -6.58749100
 H                  5.29879700   -0.92409700   -7.72565300
 H                  7.22943600   -1.96738600   -6.57099000
 H                  8.60680500   -2.20714000   -4.43825700
 H                  8.90531600   -1.61115000   -2.04669200
 H                  5.66870500    0.99665500   -2.13722100
 H                 -8.92604400    1.63095900    2.01879000
 H                 -8.63280000    3.82838500    3.13805700
 H                 -7.25644600    5.95967100    3.38896300
 H                 -5.32069400    7.31838000    2.64100900
 H                 -3.92067000    6.59306600    0.72581900
 H                 -5.68082800    2.30672700   -0.49027300
 H                 -5.90313100    4.42711000   -2.16158700
 H                 -4.58980900    3.69227800   -4.12937900
 H                 -1.03033100    3.54746100   -1.68524500
 H                 -2.34251500    4.25788000    0.26880600
 H                 -1.26625100    2.66008700   -7.77430600
 H                 -0.04275100    1.03380100   -9.24465600
 H                  1.20139800   -0.86964300   -8.18036500
 H                  1.01392500   -3.07407200   -2.44234300
 H                  2.32758300   -4.20238500   -0.69714800
 H                  5.88751300   -3.81743200   -3.10392200
 H                  4.57228000   -2.66289400   -4.85733100
 H                  3.90780500   -6.58257600   -0.78378100
 H                  5.30697100   -7.72179700    0.91759900
 H                  7.24176100   -6.56751600    1.95518400
 H                  8.61882600   -4.43400800    2.19268700
 H                  8.91122500   -2.04046100    1.60061400
 H                  5.66623000   -2.13089700   -0.99644000
 N                 -7.28285500    0.73283400   -1.89471200
 N                  0.00889700   -5.41033200   -0.59752400
 N                  7.28243800    0.29287100    2.00240100
 N                 -7.29108500   -1.88954000   -0.73069400
 N                  0.01220700   -0.60801700    5.39171100
 N                  7.28468400    2.00834200   -0.29721800
 N                 -7.29529000   -0.72651900    1.89088000
 N                  0.00226400    5.40059200    0.61019000
 N                  7.28062900   -0.28962300   -2.01310600
 N                 -7.29700600    1.89611300    0.72937500
 N                 -0.01092200    0.61306000   -5.39301300
 N                  7.28222500   -2.00676600    0.28477100
 Pd                -7.37603800    0.00368600   -0.00200000
 Pd                 7.36507300    0.00094100   -0.00539800

