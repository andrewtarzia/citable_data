%nprocshared=48
%mem=120GB
%chk=m6_l1_opt_FINAL.chk
# SP HSEH1PBE/def2svp scrf=(solvent=dmso)

m6_l1

12 1
 Pd                -7.89095200   -2.76569100    6.75891600
 Pd                -7.25495700    4.21422100   -6.73251500
 Pd                 7.89098900    2.76950000   -6.77043300
 Pd                 7.26243100   -4.20057700    6.72341700
 Pd                -0.91674800   -9.50669100   -4.94406000
 Pd                 0.92108000    9.50262400    4.94357500
 C                 -8.90224100   -0.89916500    4.70848700
 N                 -9.25361300   -1.93111900    5.49496800
 C                -10.50867600   -2.40565700    5.41993100
 C                -11.44888400   -1.86765600    4.55767900
 C                -11.10566800   -0.79473200    3.71894200
 C                -12.07147900   -0.20110800    2.77567500
 C                -13.45483300   -0.21967900    3.06689100
 C                -14.39809600    0.33156800    2.20718200
 C                -13.91947700    0.89857600    1.03399300
 O                -14.68015000    1.47892800    0.06899500
 C                -13.82741500    1.89991700   -0.90180000
 C                -14.19464000    2.54614300   -2.07414400
 C                -13.16913200    2.89102200   -2.94708100
 C                -11.81321500    2.61121700   -2.66081800
 C                -10.75786400    2.99951000   -3.61503200
 C                 -9.45826200    3.31624500   -3.18570200
 C                 -8.48802600    3.66979900   -4.10640100
 N                 -8.73906200    3.72078300   -5.42593200
 C                 -9.97522300    3.42554600   -5.86262000
 C                -10.99574800    3.06923500   -4.99750400
 C                -11.47569600    1.96125600   -1.46480100
 C                -12.49049600    1.59916900   -0.58126300
 C                -12.55166800    0.93402200    0.70476300
 C                -11.62067200    0.38433400    1.58364700
 C                 -9.78874400   -0.31592500    3.82081400
 H                 -7.87975400   -0.52728600    4.80137100
 H                -10.76268000   -3.24846800    6.06579000
 H                -12.44229700   -2.31591700    4.53104100
 H                -13.80091500   -0.65390900    4.00663000
 H                -15.46250000    0.32512200    2.44437700
 H                -15.23728400    2.77669200   -2.29545500
 H                -13.42681800    3.41859500   -3.86727400
 H                 -9.19301800    3.31591900   -2.12822000
 H                 -7.47805300    3.92723100   -3.78132700
 H                -10.14710500    3.46854400   -6.94000500
 H                -11.96917900    2.82025900   -5.42064100
 H                -10.43701900    1.71440800   -1.23743500
 H                -10.55916100    0.38668200    1.32922700
 H                 -9.44472700    0.53200800    3.22820500
 C                 -5.36985100   -4.10447900    7.50907100
 N                 -6.50954400   -3.59408000    8.00686300
 C                 -6.70933500   -3.62784800    9.33538700
 C                 -5.77875000   -4.16895000   10.20606800
 C                 -4.57629400   -4.70160300    9.71327700
 C                 -3.55369700   -5.27663500   10.60698300
 C                 -3.93372600   -5.87954400   11.82802900
 C                 -3.00337800   -6.43752700   12.69736200
 C                 -1.67053100   -6.36802400   12.31591500
 O                 -0.61822200   -6.84801600   13.02935200
 C                  0.50098100   -6.57482200   12.30857400
 C                  1.79920800   -6.89530200   12.68099400
 C                  2.81240600   -6.52104600   11.80583300
 C                  2.54514000   -5.85433600   10.58812200
 C                  3.65262000   -5.47873300    9.68953500
 C                  3.47468400   -5.39297500    8.29854800
 C                  4.53156100   -5.03168900    7.48201900
 N                  5.75038900   -4.74957000    7.97398000
 C                  5.94768800   -4.82621300    9.30112200
 C                  4.93683200   -5.18495300   10.17635400
 C                  1.22127300   -5.55048900   10.23886000
 C                  0.19256200   -5.90943500   11.10747400
 C                 -1.25012400   -5.77274900   11.11199700
 C                 -2.19852000   -5.22736000   10.24913000
 C                 -4.39492600   -4.65645500    8.32077900
 H                 -5.24193400   -4.07475800    6.42499400
 H                 -7.64084200   -3.19953800    9.71070100
 H                 -5.99413600   -4.14351300   11.27445800
 H                 -4.99113300   -5.93879800   12.09232300
 H                 -3.30511200   -6.91310400   13.63118200
 H                  2.01192300   -7.42054900   13.61269900
 H                  3.84135200   -6.77875300   12.06328400
 H                  2.51857200   -5.63153400    7.83221300
 H                  4.40450200   -4.97237800    6.39899200
 H                  6.94565600   -4.58399900    9.67149900
 H                  5.15949800   -5.20513700   11.24339500
 H                  0.99608300   -5.01487400    9.31472300
 H                 -1.88239300   -4.74482600    9.32245800
 H                 -3.50483000   -5.07384700    7.84957000
 C                 -7.78318500   -4.64767200    4.48612100
 N                 -8.08382000   -4.54933200    5.79260100
 C                 -8.49088100   -5.64925000    6.44897600
 C                 -8.60569500   -6.87899700    5.82361900
 C                 -8.28842500   -7.01155600    4.46199200
 C                 -8.38204900   -8.30697300    3.76334800
 C                 -9.32792400   -9.27221000    4.17863000
 C                 -9.44893100  -10.50592300    3.54946600
 C                 -8.58601300  -10.75570300    2.49126900
 O                 -8.54390200  -11.89999700    1.75947200
 C                 -7.56218100  -11.73730800    0.83406600
 C                 -7.19153200  -12.67259400   -0.12239500
 C                 -6.16510400  -12.31101500   -0.98745200
 C                 -5.51660500  -11.05774900   -0.90084200
 C                 -4.43296600  -10.71599100   -1.84084800
 C                 -3.39251900   -9.84406500   -1.47904200
 C                 -2.39742500   -9.53243400   -2.38838300
 N                 -2.38468700  -10.03669900   -3.63437800
 C                 -3.36533300  -10.87834000   -4.00356400
 C                 -4.38925100  -11.23769800   -3.14400400
 C                 -5.91167300  -10.13746500    0.08087100
 C                 -6.94494800  -10.47801300    0.95145200
 C                 -7.62862700   -9.82380700    2.04907100
 C                 -7.52914500   -8.58852800    2.68627500
 C                 -7.87096100   -5.84545500    3.79904500
 H                 -7.46892400   -3.73435400    3.97679500
 H                 -8.72214400   -5.53980000    7.51039300
 H                 -8.91885100   -7.73585500    6.42050300
 H                -10.00710000   -9.04052400    5.00118600
 H                -10.18891600  -11.24119000    3.86734100
 H                 -7.67915800  -13.64589700   -0.18706300
 H                 -5.83922200  -13.03463500   -1.73696200
 H                 -3.33171800   -9.41680400   -0.47783300
 H                 -1.57944100   -8.86348200   -2.11342600
 H                 -3.33144300  -11.26730800   -5.02316700
 H                 -5.16573900  -11.90669200   -3.51552100
 H                 -5.43987000   -9.15537800    0.14785300
 H                 -6.77598900   -7.86579300    2.36675800
 H                 -7.63319000   -5.85060200    2.73506300
 C                 -6.48625500   -0.36383000    7.74743000
 N                 -7.68199900   -0.97691200    7.71253900
 C                 -8.73323700   -0.38139800    8.30110900
 C                 -8.62406300    0.84290600    8.93815300
 C                 -7.38841400    1.50894100    8.98136200
 C                 -7.23219500    2.82035700    9.63761400
 C                 -8.05496600    3.17059000   10.73263900
 C                 -7.93849000    4.39408000   11.38235000
 C                 -6.97645500    5.27010400   10.89865500
 O                 -6.71599300    6.50907000   11.39257700
 C                 -5.71279400    7.02890300   10.63753800
 C                 -5.14323400    8.28362600   10.80520600
 C                 -4.12923300    8.63454300    9.92132900
 C                 -3.68482600    7.76068300    8.90266900
 C                 -2.60326400    8.17457500    7.98953000
 C                 -1.73653800    7.23704200    7.40322000
 C                 -0.73508800    7.65576200    6.54542300
 N                 -0.55317400    8.95165900    6.23763200
 C                 -1.36734900    9.86903600    6.78687000
 C                 -2.38651200    9.52184200    7.65719700
 C                 -4.27797900    6.49739400    8.76395700
 C                 -5.30187700    6.13076800    9.63499400
 C                 -6.14094300    4.96210700    9.80886600
 C                 -6.26624400    3.72665500    9.17678600
 C                 -6.30561400    0.86216400    8.36263000
 H                 -5.64858600   -0.87836900    7.27232800
 H                 -9.69270200   -0.89988200    8.25066000
 H                 -9.52192800    1.28105200    9.37447000
 H                 -8.78886300    2.45342000   11.10489900
 H                 -8.56751500    4.65215500   12.23497100
 H                 -5.47230800    8.95806700   11.59644100
 H                 -3.64906700    9.60723400   10.04304400
 H                 -1.81502800    6.17306300    7.62727300
 H                 -0.05028100    6.93527200    6.09336600
 H                 -1.20062000   10.91257300    6.51229700
 H                 -3.02435500   10.31440300    8.04896000
 H                 -3.96578500    5.81814800    7.96853500
 H                 -5.63818000    3.48419000    8.31755100
 H                 -5.30421000    1.29304400    8.37404500
 C                 -4.63195400    3.96543700   -8.06414700
 N                 -5.75649900    4.70082700   -8.02428400
 C                 -5.86228400    5.76871000   -8.83332400
 C                 -4.84905500    6.13643100   -9.70216700
 C                 -3.66004100    5.39099500   -9.75730100
 C                 -2.55225800    5.75715900  -10.65940900
 C                 -2.81932100    6.41469900  -11.88215500
 C                 -1.80651100    6.77697400  -12.76268800
 C                 -0.50740800    6.47027400  -12.38193500
 O                  0.61253800    6.74781500  -13.09976100
 C                  1.66625100    6.28882700  -12.37475600
 C                  2.99994900    6.36651700  -12.75153000
 C                  3.93176300    5.82808100  -11.87160500
 C                  3.55213000    5.23590700  -10.64521800
 C                  4.57625600    4.67635600   -9.74351900
 C                  4.39384400    4.64472600   -8.35094600
 C                  5.37002000    4.10359200   -7.53339600
 N                  6.51110900    3.59071200   -8.02512400
 C                  6.71226700    3.61303200   -9.35368000
 C                  5.78119600    4.14352600  -10.23014800
 C                  2.19611700    5.17900700  -10.29170600
 C                  1.24621400    5.70489000  -11.16512400
 C                 -0.19842100    5.82039600  -11.17247600
 C                 -1.22831800    5.45715100  -10.30692500
 C                 -3.57826200    4.27744100   -8.90474500
 H                 -4.58164300    3.09633300   -7.40516600
 H                 -6.78431900    6.35035400   -8.77397700
 H                 -4.99057500    7.02761100  -10.31396500
 H                 -3.85231400    6.62465500  -12.16553400
 H                 -2.02101100    7.27301200  -13.70986100
 H                  3.30100900    6.83286500  -13.69020700
 H                  4.98969300    5.89357200  -12.13234900
 H                  3.50159600    5.06333600   -7.88490000
 H                  5.24182700    4.08386700   -6.44942000
 H                  7.64499600    3.18329100   -9.72421800
 H                  5.99775100    4.10871800  -11.29806500
 H                  1.88083900    4.70412200   -9.36075700
 H                 -0.99944900    4.97253500   -9.35594600
 H                 -2.70447900    3.62548400   -8.90384600
 C                  8.90786400    0.91006800   -4.71396200
 N                  9.25433000    1.94257900   -5.50193100
 C                 10.50574200    2.42643800   -5.42443200
 C                 11.44706800    1.89727000   -4.55801800
 C                 11.10878300    0.82427300   -3.71751200
 C                 12.07643300    0.24206800   -2.76908000
 C                 13.45954000    0.26266300   -3.06111700
 C                 14.40507000   -0.27503300   -2.19537200
 C                 13.92874500   -0.83036600   -1.01570100
 O                 14.69221600   -1.39493900   -0.04354600
 C                 13.84181000   -1.80533700    0.93367800
 C                 14.21172400   -2.44049800    2.11122300
 C                 13.18699100   -2.78938900    2.98341300
 C                 11.83084400   -2.50800300    2.69983300
 C                 10.77679900   -2.89599000    3.65549900
 C                  9.57874000   -2.17067200    3.76457900
 C                  8.60453500   -2.56713800    4.66333500
 N                  8.75650700   -3.64336300    5.45433500
 C                  9.89606200   -4.35132500    5.37380000
 C                 10.91408800   -4.00991400    4.49982300
 C                 11.49149300   -1.86011800    1.50319700
 C                 12.50386700   -1.51046300    0.61186600
 C                 12.56116700   -0.86751600   -0.68562400
 C                 11.62775300   -0.33181000   -1.57072000
 C                  9.79581600    0.33556600   -3.82183600
 H                  7.88829700    0.53076700   -4.80899800
 H                 10.75490600    3.26983600   -6.07132000
 H                 12.43732600    2.35228000   -4.52886900
 H                 13.80333100    0.68790200   -4.00582800
 H                 15.46948200   -0.26675300   -2.43248300
 H                 15.25679300   -2.64966800    2.34212000
 H                 13.45040200   -3.27123500    3.92674900
 H                  9.40025900   -1.27488600    3.16935500
 H                  7.67420800   -2.00405800    4.76078700
 H                  9.98900700   -5.22318100    6.02446100
 H                 11.80207800   -4.64149500    4.46779800
 H                 10.44825800   -1.65612600    1.25502300
 H                 10.56619300   -0.33484500   -1.31633000
 H                  9.45655300   -0.51322000   -3.22761800
 C                  6.49102100    0.36029800   -7.74589900
 N                  7.68514500    0.97695800   -7.71747700
 C                  8.73575000    0.38309200   -8.30902600
 C                  8.62707600   -0.84241200   -8.94369000
 C                  7.39293900   -1.51172800   -8.98117000
 C                  7.23778400   -2.82402900   -9.63589400
 C                  8.06300200   -3.17565400  -10.72866700
 C                  7.94956300   -4.40076200  -11.37572300
 C                  6.98747500   -5.27676800  -10.89217000
 O                  6.73065900   -6.51759000  -11.38307800
 C                  5.72602400   -7.03680900  -10.62954500
 C                  5.16039300   -8.29362500  -10.79449400
 C                  4.14420900   -8.64375000   -9.91287700
 C                  3.69329100   -7.76658200   -8.90002800
 C                  2.60873000   -8.17990400   -7.99032400
 C                  1.72812500   -7.24464700   -7.42193200
 C                  0.72489800   -7.66385400   -6.56638100
 N                  0.55376700   -8.95775100   -6.24432200
 C                  1.38064800   -9.87314600   -6.77795200
 C                  2.40274000   -9.52540600   -7.64458100
 C                  4.28266100   -6.50130200   -8.76362700
 C                  5.30979400   -6.13598900   -9.63155900
 C                  6.14821200   -4.96684300   -9.80578300
 C                  6.27120000   -3.73008500   -9.17578900
 C                  6.31115000   -0.86716800   -8.35833200
 H                  5.65411900    0.87307900   -7.26750400
 H                  9.69416900    0.90434900   -8.26439500
 H                  9.52429200   -1.27848300   -9.38326100
 H                  8.79698300   -2.45860900  -11.10094500
 H                  8.58102200   -4.66020800  -12.22612500
 H                  5.49415200   -8.97025000  -11.58189200
 H                  3.66699200   -9.61821700  -10.03214400
 H                  1.79731400   -6.18274800   -7.65849300
 H                  0.02965600   -6.94534600   -6.12798400
 H                  1.22195500  -10.91519100   -6.49321500
 H                  3.05210600  -10.31521600   -8.02264900
 H                  3.96509400   -5.81984100   -7.97216700
 H                  5.64092000   -3.48659600   -8.31849100
 H                  5.31114100   -1.30121000   -8.36468900
 C                  7.76721800    4.65310600   -4.50058600
 N                  8.07612000    4.55453600   -5.80527800
 C                  8.48882600    5.65400400   -6.45887900
 C                  8.60388000    6.88288200   -5.83162200
 C                  8.28076700    7.01513700   -4.47129300
 C                  8.37950800    8.30806600   -3.76878100
 C                  9.33145700    9.26954400   -4.17894700
 C                  9.46116400   10.49859500   -3.54242000
 C                  8.60095700   10.74761700   -2.48183900
 O                  8.56883600   11.88656200   -1.74131000
 C                  7.58711700   11.72450400   -0.81573500
 C                  7.22620500   12.65466600    0.14936000
 C                  6.19794700   12.29462900    1.01289500
 C                  5.53736600   11.04848600    0.91570800
 C                  4.45062600   10.70927800    1.85316900
 C                  3.39840600    9.85582600    1.48206200
 C                  2.40036000    9.54652200    2.38886100
 N                  2.39498000   10.03581800    3.64088500
 C                  3.38732300   10.85941100    4.01931900
 C                  4.41512500   11.21532700    3.16272300
 C                  5.92304000   10.13316800   -0.07437300
 C                  6.95886800   10.47162100   -0.94275300
 C                  7.63673200    9.82001700   -2.04549700
 C                  7.52824900    8.58954700   -2.69035800
 C                  7.85465600    5.84998300   -3.81214700
 H                  7.44718700    3.74069000   -3.99306500
 H                  8.72558900    5.54472900   -7.51915200
 H                  8.92287800    7.73921600   -6.42605800
 H                 10.00930300    9.03805200   -5.00262700
 H                 10.20645500   11.23040500   -3.85584200
 H                  7.72271500   13.62291200    0.22188800
 H                  5.87988300   13.01450100    1.76931100
 H                  3.33089600    9.44201700    0.47574200
 H                  1.57352800    8.89200100    2.10662000
 H                  3.35956500   11.23595900    5.04372000
 H                  5.20146000   11.86865500    3.54131800
 H                  5.44255000    9.15586600   -0.14898200
 H                  6.77058500    7.87001000   -2.37436700
 H                  7.61020900    5.85473100   -2.74966200
 C                 -7.23700500    1.26825700   -6.92996000
 N                 -7.44510700    2.39565400   -7.63194400
 C                 -7.77000700    2.29895000   -8.93250100
 C                 -7.89562900    1.07681900   -9.57080700
 C                 -7.67661400   -0.11602900   -8.86253100
 C                 -7.78865900   -1.43779200   -9.50711600
 C                 -8.67059900   -1.62608000  -10.59592400
 C                 -8.80939700   -2.85537900  -11.23007700
 C                 -8.02982700   -3.89858900  -10.74970000
 O                 -8.01798600   -5.16642900  -11.23880300
 C                 -7.12208500   -5.86545500  -10.49359900
 C                 -6.80150000   -7.20531600  -10.66370700
 C                 -5.85022300   -7.73642300   -9.80023600
 C                 -5.23915800   -6.96175700   -8.78767600
 C                 -4.23750800   -7.56826700   -7.89124100
 C                 -4.05427400   -7.10979000   -6.57596500
 C                 -3.10009700   -7.69307000   -5.76147800
 N                 -2.31992300   -8.70293400   -6.18462000
 C                 -2.47756800   -9.16256700   -7.43772200
 C                 -3.41495000   -8.62853600   -8.30533900
 C                 -5.58940800   -5.61187800   -8.63948500
 C                 -6.53429000   -5.05889700   -9.50135500
 C                 -7.13744700   -3.75224000   -9.67140900
 C                 -7.01923900   -2.51441500   -9.04258600
 C                 -7.34111700    0.01407500   -7.50459600
 H                 -6.98630400    1.37765700   -5.87279100
 H                 -7.92617400    3.23156900   -9.47852900
 H                 -8.14020600    1.07099900  -10.63319700
 H                 -9.28687100   -0.79311700  -10.93915400
 H                 -9.50144000   -2.99660600  -12.06095300
 H                 -7.27818800   -7.81220900  -11.43422100
 H                 -5.59399900   -8.79298700   -9.89759200
 H                 -4.66860400   -6.30966800   -6.16246600
 H                 -2.95754800   -7.34765100   -4.73532600
 H                 -1.82517200   -9.97836300   -7.75529600
 H                 -3.47411400   -9.03339000   -9.31584000
 H                 -5.11090600   -4.99105300   -7.87982100
 H                 -6.31520700   -2.38341800   -8.21872300
 H                 -7.18306800   -0.85963900   -6.87211400
 C                 -5.88624500    6.35820400   -5.23347600
 N                 -7.05071100    6.02256200   -5.81531300
 C                 -8.07674700    6.88879800   -5.75934300
 C                 -7.97324800    8.11154600   -5.11827700
 C                 -6.77155200    8.48227800   -4.49299800
 C                 -6.62613300    9.76675800   -3.78312500
 C                 -7.37979500   10.89281400   -4.18663400
 C                 -7.27527300   12.12076400   -3.54338400
 C                 -6.39654100   12.19244100   -2.47118600
 O                 -6.16292700   13.29587100   -1.71319500
 C                 -5.24725200   12.93957700   -0.77437500
 C                 -4.72599500   13.77600600    0.20311100
 C                 -3.78904500   13.22147000    1.06749200
 C                 -3.38549100   11.86988400    0.97139600
 C                 -2.38764000   11.32469900    1.91065000
 C                 -2.38393200    9.96818700    2.27598800
 C                 -1.42806200    9.48568300    3.15224300
 N                 -0.47850600   10.27686300    3.68093700
 C                 -0.46307100   11.57987400    3.35191000
 C                 -1.39043400   12.13053200    2.48386300
 C                 -3.93604800   11.05006800   -0.02447600
 C                 -4.87018000   11.59014400   -0.90613500
 C                 -5.63016800   11.09512300   -2.03636400
 C                 -5.74232300    9.87489200   -2.69958600
 C                 -5.71379900    7.56116100   -4.57228100
 H                 -5.06636000    5.64069200   -5.30671200
 H                 -9.01172500    6.58624000   -6.23523500
 H                 -8.85192800    8.75611800   -5.08909900
 H                 -8.04827600   10.81290700   -5.04582200
 H                 -7.85122000   12.98815700   -3.86765500
 H                 -5.04134200   14.81627200    0.29078200
 H                 -3.37757100   13.85057200    1.85892600
 H                 -3.13959900    9.27683700    1.90250000
 H                 -1.42367700    8.43467800    3.44828800
 H                  0.32272700   12.19616200    3.79324100
 H                 -1.30639900   13.19019900    2.24181700
 H                 -3.61901200   10.01098100   -0.13101400
 H                 -5.17027400    9.01070300   -2.35694700
 H                 -4.73640300    7.77910600   -4.14134100
 C                  6.79309000   -6.02544100    4.45277000
 N                  7.11046400   -5.98895200    5.75842300
 C                  7.30537700   -7.14762600    6.41085100
 C                  7.18710200   -8.37521800    5.78167000
 C                  6.84889300   -8.44181900    4.42009800
 C                  6.70104200   -9.72888900    3.71539100
 C                  7.44663000  -10.85651900    4.12933500
 C                  7.34000500  -12.08727100    3.49180500
 C                  6.45364400  -12.16509200    2.42635200
 O                  6.20435400  -13.27652800    1.68522500
 C                  5.27861000  -12.92683100    0.75378500
 C                  4.74080500  -13.77224500   -0.20688500
 C                  3.79627400  -13.22313900   -1.06654900
 C                  3.40288300  -11.86778800   -0.98310300
 C                  2.40033000  -11.32679800   -1.91967100
 C                  2.41117500   -9.97723900   -2.30964800
 C                  1.45236400   -9.49768000   -3.18417900
 N                  0.48580300  -10.28495700   -3.68730300
 C                  0.45505700  -11.58102900   -3.33279100
 C                  1.38505700  -12.12888500   -2.46571700
 C                  3.96968400  -11.03895700   -0.00398900
 C                  4.91031500  -11.57387500    0.87384400
 C                  5.69062600  -11.06794000    1.98520100
 C                  5.81838000   -9.83998900    2.63118200
 C                  6.65472700   -7.21601300    3.76193700
 H                  6.65334100   -5.06778300    3.94754600
 H                  7.55542500   -7.08682000    7.47195800
 H                  7.33575500   -9.27725600    6.37552100
 H                  8.15109600  -10.76097000    4.95762000
 H                  7.92769600  -12.94943400    3.80887700
 H                  5.04937800  -14.81526100   -0.28556400
 H                  3.37104400  -13.85949800   -1.84483100
 H                  3.18084000   -9.29004100   -1.95766400
 H                  1.45907200   -8.45263100   -3.50032800
 H                 -0.34572000  -12.19377600   -3.75136500
 H                  1.28867700  -13.18226500   -2.20188900
 H                  3.66032900   -9.99654400    0.09256500
 H                  5.21555100   -8.98778800    2.31201700
 H                  6.41845000   -7.17296400    2.69857000
 C                  6.34880500   -1.56544400    7.69496400
 N                  7.39926100   -2.40394000    7.67515200
 C                  8.53977400   -2.02896400    8.27947000
 C                  8.66671600   -0.80705800    8.91792300
 C                  7.58740900    0.09115000    8.94622400
 C                  7.68657900    1.40825200    9.60234800
 C                  8.55367900    1.59094700   10.70400300
 C                  8.67786700    2.81499200   11.35114900
 C                  7.91255300    3.86328200   10.85903700
 O                  7.90104300    5.13085000   11.34864600
 C                  7.02748000    5.83799500   10.58472500
 C                  6.71278900    7.17973900   10.75091400
 C                  5.78368000    7.71922900    9.86868600
 C                  5.18858200    6.95057800    8.84223000
 C                  4.20586100    7.56427600    7.92981900
 C                  4.04827200    7.11508400    6.60813700
 C                  3.10754900    7.70215600    5.78072200
 N                  2.31670300    8.70624600    6.19724000
 C                  2.45155400    9.15844400    7.45566200
 C                  3.37495800    8.62104900    8.33609100
 C                  5.53296500    5.59880100    8.69784600
 C                  6.45558000    5.03746400    9.57839000
 C                  7.04007600    3.72432100    9.76357500
 C                  6.92289200    2.48703300    9.13339400
 C                  6.40668500   -0.32787400    8.31075600
 H                  5.43235300   -1.90252600    7.20635100
 H                  9.37807100   -2.72749400    8.24081600
 H                  9.62755400   -0.55613500    9.36784400
 H                  9.12827800    0.74367200   11.08272900
 H                  9.34013300    2.94570700   12.20762900
 H                  7.17680800    7.78150200   11.53310500
 H                  5.53185000    8.77716500    9.96301200
 H                  4.67213600    6.31935700    6.20036800
 H                  2.98431300    7.36514700    4.74950500
 H                  1.79131100    9.97021700    7.76701000
 H                  3.41585700    9.01909500    9.35021100
 H                  5.06542900    4.98302100    7.92725700
 H                  6.26686200    2.37246500    8.26847800
 H                  5.51089900    0.29347100    8.30865700

