%nprocshared=48
%mem=120GB
%oldchk=cis_l1_lb_3.chk
%chk=cis_l1_lb_4.chk
# OPT=restart HSEH1PBE/def2svp scrf=(solvent=dmso)

cis_l1_lb

4 1
 Pd                 7.22965200   -1.56244400    0.36151300
 Pd                -7.22951700   -1.56238200    0.36009100
 C                  4.86509200   -2.13111400   -1.30887400
 N                  5.97186500   -2.68624400   -0.79232700
 C                  6.18477600   -3.99897700   -0.97105700
 C                  5.28113500   -4.80680000   -1.64385400
 C                  4.09481500   -4.25793300   -2.15449200
 C                  3.05409400   -5.08646600   -2.79251200
 C                  3.38892100   -6.23051800   -3.55134100
 C                  2.41753200   -7.03826400   -4.13570700
 C                  1.09202300   -6.67937800   -3.92614700
 O                  0.00018800   -7.35041600   -4.38337100
 C                 -1.09191100   -6.67925200   -3.92694800
 C                 -2.41730700   -7.03804200   -4.13739600
 C                 -3.38903300   -6.23019100   -3.55372700
 C                 -3.05461500   -5.08601000   -2.79492300
 C                 -4.09559900   -4.25745000   -2.15736200
 C                 -3.91922500   -2.87425100   -1.99343000
 C                 -4.86611800   -2.13062600   -1.31205900
 N                 -5.97256700   -2.68591400   -0.79497000
 C                 -6.18533100   -3.99869500   -0.97349600
 C                 -5.28173100   -4.80645100   -1.64644000
 C                 -1.70858200   -4.74345400   -2.61227600
 C                 -0.72317800   -5.54998600   -3.17433400
 C                  0.72287500   -5.55012400   -3.17371100
 C                  1.70795400   -4.74391000   -2.61062400
 C                  3.91816200   -2.87481700   -1.99011800
 H                  4.73349700   -1.05668200   -1.16612300
 H                  7.09862300   -4.41340500   -0.54056200
 H                  5.49291500   -5.87358800   -1.72556900
 H                  4.44006300   -6.48246600   -3.70650100
 H                  2.68334100   -7.91411500   -4.72877700
 H                 -2.68277800   -7.91393200   -4.73056000
 H                 -4.44009600   -6.48212800   -3.70944800
 H                 -3.04624800   -2.36293800   -2.40012300
 H                 -4.73477600   -1.05612600   -1.16953100
 H                 -7.09896400   -4.41322100   -0.54263800
 H                 -5.49330000   -5.87330000   -1.72792200
 H                 -1.43119600   -3.87950300   -2.00567700
 H                  1.43020700   -3.88027400   -2.00373600
 H                  3.04498700   -2.36357300   -2.39646900
 C                  4.88476600   -1.14344400    2.10283400
 N                  5.98133000   -1.89974400    1.94303400
 C                  6.18295200   -2.93730900    2.76956400
 C                  5.27724600   -3.27386100    3.76344300
 C                  4.09972400   -2.52678500    3.92056600
 C                  3.05565800   -2.90233800    4.89301600
 C                  3.38893200   -3.49192400    6.13313100
 C                  2.41639300   -3.88164000    7.04919000
 C                  1.09168600   -3.68074600    6.68292800
 O                 -0.00025500   -4.02058100    7.41986700
 C                 -1.09187100   -3.68114600    6.68226900
 C                 -2.41672100   -3.88245500    7.04777600
 C                 -3.38885100   -3.49305100    6.13115800
 C                 -3.05507000   -2.90350600    4.89115400
 C                 -4.09880900   -2.52809300    3.91830100
 C                 -3.93417500   -1.42276100    3.06823200
 C                 -4.88352100   -1.14475600    2.10039900
 N                 -5.98037400   -1.90067600    1.94080000
 C                 -6.18216400   -2.93819700    2.76735500
 C                 -5.27648000   -3.27494000    3.76118000
 C                 -1.70942300   -2.70102500    4.55778400
 C                 -0.72305600   -3.10060900    5.45544200
 C                  0.72338500   -3.10029100    5.45591000
 C                  1.71014500   -2.70016500    4.55893600
 C                  3.93540100   -1.42125700    3.07069600
 H                  4.76312100   -0.29242800    1.42908900
 H                  7.08805300   -3.52637700    2.60878000
 H                  5.47880000   -4.15049900    4.38034900
 H                  4.43985100   -3.62950300    6.39594900
 H                  2.68037200   -4.32721000    8.00902900
 H                 -2.68111200   -4.32805600    8.00748700
 H                 -4.43987100   -3.63085900    6.39344400
 H                 -3.06772100   -0.76641400    3.15566000
 H                 -4.76164500   -0.29384000    1.42657300
 H                 -7.08744100   -3.52702900    2.60671000
 H                 -5.47825500   -4.15148600    4.37814300
 H                 -1.43362100   -2.27351100    3.59203700
 H                  1.43476700   -2.27243200    3.59317300
 H                  3.06911700   -0.76470100    3.15824200
 C                 10.42269000   -0.19369700    2.84910600
 C                  9.57180600   -0.94966400    2.08793900
 N                  8.44275800   -0.41889400    1.54159900
 C                  8.11835900    0.83922700    1.78550400
 C                  8.89534000    1.68946400    2.60295000
 C                  8.51856800    3.04402800    2.88214000
 C                  7.20543700    3.60337500    2.48287100
 C                  7.13046800    4.83679600    1.81586200
 C                  5.90358900    5.41232100    1.51175900
 C                  4.70319000    4.77010200    1.87229700
 C                  3.44557000    5.38410100    1.60547800
 C                  2.37467100    5.93103300    1.40470800
 C                  1.15005000    6.64230900    1.18706900
 N                  0.00067600    5.97946600    1.36889000
 C                 -1.14904100    6.64195300    1.18794100
 C                 -2.37328600    5.93030500    1.40646300
 C                 -3.44389900    5.38306800    1.60793100
 C                 -4.70123900    4.76879200    1.87542800
 C                 -5.90197300    5.41146900    1.51682000
 C                 -7.12856600    4.83564300    1.82148900
 C                 -7.20291800    3.60145500    2.48714800
 C                 -8.51573200    3.04172600    2.88692800
 C                 -9.40871700    3.82962900    3.59290000
 C                -10.63495700    3.31797000    4.07100200
 C                -10.98297400    2.00354300    3.86141000
 C                -10.11945500    1.16285000    3.12017000
 C                -10.42064200   -0.19548400    2.85100600
 C                 -9.57052600   -0.95062000    2.08816700
 N                 -8.44170400   -0.41943900    1.54174300
 C                 -8.11680700    0.83827800    1.78708400
 C                 -8.89299900    1.68761300    2.60623500
 C                 -6.00262400    2.95684200    2.83332700
 C                 -4.77191900    3.52651900    2.53359000
 C                 -1.19769000    7.99845000    0.81455900
 C                 -0.00002000    8.67847400    0.62914500
 C                  1.19800000    7.99881700    0.81364500
 C                  4.77447800    3.52856400    2.53178300
 C                  6.00545600    2.95919700    2.83096200
 C                  9.41232900    3.83270100    3.58626800
 C                 10.63885700    3.32141400    4.06403500
 C                 10.98639000    2.00661000    3.85600000
 C                 10.12207600    1.16509600    3.11662000
 H                 11.33440300   -0.64225100    3.24702300
 H                  9.76633300   -2.00308500    1.88139000
 H                  7.21663100    1.22036300    1.30303700
 H                  8.05060200    5.35005000    1.52731300
 H                  5.86497400    6.37323000    0.99528500
 H                 -5.86383600    6.37295900    1.00139300
 H                 -8.04897200    5.34924700    1.53443500
 H                 -9.14106500    4.86378700    3.82038500
 H                -11.29993100    3.97584200    4.63474300
 H                -11.91796000    1.59865000    4.25377600
 H                -11.33217400   -0.64434600    3.24899000
 H                 -9.76550000   -2.00369900    1.88030100
 H                 -7.21538500    1.21992000    1.30445500
 H                 -6.03036700    2.01248700    3.38212300
 H                 -3.85085300    3.02069100    2.82914300
 H                 -2.16034200    8.49407800    0.68250200
 H                 -0.00029000    9.73256800    0.34367500
 H                  2.16039900    8.49474000    0.68085300
 H                  3.85368700    3.02310800    2.82882200
 H                  6.03370800    2.01548400    3.38084000
 H                  9.14509200    4.86721800    3.81260100
 H                 11.30445100    3.97990900    4.62631600
 H                 11.92158800    1.60203800    4.24819200
 C                 10.37905700   -1.70351300   -2.53081200
 C                  9.59174700   -1.89394100   -1.42625900
 N                  8.42088600   -1.21816800   -1.26117300
 C                  8.02699300   -0.34536100   -2.17255400
 C                  8.76756700   -0.06614700   -3.34200500
 C                  8.30822600    0.85645300   -4.33753900
 C                  7.05170100    1.62654200   -4.16744500
 C                  6.02189100    1.50581700   -5.11289100
 C                  4.83334400    2.21056500   -4.96849900
 C                  4.64528000    3.07460200   -3.87420700
 C                  3.42009800    3.78629700   -3.71065900
 C                  2.36832200    4.38314600   -3.56105400
 C                  1.14899500    5.11106400   -3.35942600
 N                 -0.00088100    4.45637600   -3.56306700
 C                 -1.15028100    5.11182900   -3.35919400
 C                 -2.37011600    4.38469800   -3.56059600
 C                 -3.42226100    3.78845700   -3.71002200
 C                 -4.64778500    3.07731800   -3.87344700
 C                 -4.83658700    2.21398500   -4.96816500
 C                 -6.02537400    1.50961200   -5.11239900
 C                 -7.05469300    1.63003000   -4.16638300
 C                 -8.31134400    0.86008400   -4.33626000
 C                 -9.06301300    1.00062800   -5.48860000
 C                -10.26924200    0.28952800   -5.67805200
 C                -10.73544500   -0.58231200   -4.72140600
 C                 -9.98764400   -0.78518700   -3.53722000
 C                -10.38029000   -1.70223800   -2.53071100
 C                 -9.59231300   -1.89339600   -1.42675800
 N                 -8.42174300   -1.21714600   -1.26157700
 C                 -8.02871200   -0.34328000   -2.17231100
 C                 -8.76990100   -0.06343800   -3.34122000
 C                 -6.86972900    2.50231000   -3.08142800
 C                 -5.68891100    3.22023200   -2.93705500
 C                 -1.19863900    6.45806000   -2.95440200
 C                  0.00008800    7.13338500   -2.75705700
 C                  1.19832700    6.45726500   -2.95464400
 C                  5.68686700    3.21775500   -2.93836000
 C                  6.86743200    2.49945400   -3.08287400
 C                  9.05919300    0.99619300   -5.49044000
 C                 10.26549300    0.28524100   -5.67997500
 C                 10.73254300   -0.58555700   -4.72279100
 C                  9.98548300   -0.78760200   -3.53799700
 H                 11.30458000   -2.27153300   -2.63801100
 H                  9.86248400   -2.60537300   -0.64524400
 H                  7.07494700    0.15826200   -1.99855900
 H                  6.14840100    0.83149700   -5.96299900
 H                  4.03568700    2.09241100   -5.70424800
 H                 -4.03930000    2.09605100   -5.70435000
 H                 -6.15245100    0.83580200   -5.96282700
 H                 -8.72699300    1.69518700   -6.26169600
 H                -10.83716500    0.44508100   -6.59790800
 H                -11.66927000   -1.12950800   -4.86520100
 H                -11.30561000   -2.27057300   -2.63798400
 H                 -9.86228700   -2.60575200   -0.64632200
 H                 -7.07693200    0.16077800   -1.99811200
 H                 -7.66915000    2.63621700   -2.34851800
 H                 -5.56324300    3.89886700   -2.09108200
 H                 -2.16057600    6.94985300   -2.80389000
 H                  0.00046700    8.18044300   -2.44818800
 H                  2.16062000    6.94841600   -2.80431700
 H                  5.56177000    3.89691000   -2.09271900
 H                  7.66722900    2.63359900   -2.35041600
 H                  8.72259000    1.69005800   -6.26390600
 H                 10.83283300    0.44013100   -6.60030200
 H                 11.66648100   -1.13255100   -4.86660900

