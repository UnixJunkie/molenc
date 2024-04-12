
(* the first atomic number (0) is FAKE but necessary for tabulation *)
let anums = A.of_list (L.range 0 `To 118)

(* associate a prime number to each atomic number,
   1st elt. also for tabulation reasons *)
let primes =
  [|-1;2;3;5;7;11;13;17;19;23;29;31;37;41;43;47;53;59;61;67;71;73;
    79;83;89;97;101;103;107;109;113;127;131;137;139;149;151;157;163;
    167;173;179;181;191;193;197;199;211;223;227;229;233;239;241;251;
    257;263;269;271;277;281;283;293;307;311;313;317;331;337;347;349;
    353;359;367;373;379;383;389;397;401;409;419;421;431;433;439;443;
    449;457;461;463;467;479;487;491;499;503;509;521;523;541;547;557;
    563;569;571;577;587;593;599;601;607;613;617;619;631;641;643;647|]

(* chemical symbols; 1st elt. is also for tabulation reasons only *)
let symbols =
  [|"";"H";"He";"Li";"Be";"B";"C";"N";"O";"F";"Ne";"Na";"Mg";"Al";"Si";"P";
    "S";"Cl";"Ar";"K";"Ca";"Sc";"Ti";"V";"Cr";"Mn";"Fe";"Co";"Ni";"Cu";"Zn";
    "Ga";"Ge";"As";"Se";"Br";"Kr";"Rb";"Sr";"Y";"Zr";"Nb";"Mo";"Tc";"Ru";"Rh";
    "Pd";"Ag";"Cd";"In";"Sn";"Sb";"Te";"I";"Xe";"Cs";"Ba";"La";"Ce";"Pr";"Nd";
    "Pm";"Sm";"Eu";"Gd";"Tb";"Dy";"Ho";"Er";"Tm";"Yb";"Lu";"Hf";"Ta";"W";"Re";
    "Os";"Ir";"Pt";"Au";"Hg";"Tl";"Pb";"Bi";"Po";"At";"Rn";"Fr";"Ra";"Ac";"Th";
    "Pa";"U";"Np";"Pu";"Am";"Cm";"Bk";"Cf";"Es";"Fm";"Md";"No";"Lr";"Rf";"Db";
    "Sg";"Bh";"Hs";"Mt";"Ds";"Rg";"Cn";"Nh";"Fl";"Mc";"Lv";"Ts";"Og"|]
