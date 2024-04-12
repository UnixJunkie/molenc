
module A = BatArray
module Ht = BatHashtbl
module L = BatList
module Log = Dolog.Log

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

let elements_regexp = Str.regexp "He|Li|Be|Ne|Na|Mg|Al|Si|Cl|Ar|Ca|Sc|Ti|Cr|Mn|Fe|Co|Ni|Cu|Zn|Ga|Ge|As|Se|Br|Kr|Rb|Sr|Zr|Nb|Mo|Tc|Ru|Rh|Pd|Ag|Cd|In|Sn|Sb|Te|Xe|Cs|Ba|La|Ce|Pr|Nd|Pm|Sm|Eu|Gd|Tb|Dy|Ho|Er|Tm|Yb|Lu|Hf|Ta|Re|Os|Ir|Pt|Au|Hg|Tl|Pb|Bi|Po|At|Rn|Fr|Ra|Ac|Th|Pa|Np|Pu|Am|Cm|Bk|Cf|Es|Fm|Md|No|Lr|Rf|Db|Sg|Bh|Hs|Mt|Ds|Rg|Cn|Nh|Fl|Mc|Lv|Ts|Og|H|B|C|N|O|F|P|S|K|V|Y|I|W|U"

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

let symbol2anum =
  Ht.of_list [("H",1);
              ("He",2);
              ("Li",3);
              ("Be",4);
              ("B",5);
              ("C",6);
              ("N",7);
              ("O",8);
              ("F",9);
              ("Ne",10);
              ("Na",11);
              ("Mg",12);
              ("Al",13);
              ("Si",14);
              ("P",15);
              ("S",16);
              ("Cl",17);
              ("Ar",18);
              ("K",19);
              ("Ca",20);
              ("Sc",21);
              ("Ti",22);
              ("V",23);
              ("Cr",24);
              ("Mn",25);
              ("Fe",26);
              ("Co",27);
              ("Ni",28);
              ("Cu",29);
              ("Zn",30);
              ("Ga",31);
              ("Ge",32);
              ("As",33);
              ("Se",34);
              ("Br",35);
              ("Kr",36);
              ("Rb",37);
              ("Sr",38);
              ("Y",39);
              ("Zr",40);
              ("Nb",41);
              ("Mo",42);
              ("Tc",43);
              ("Ru",44);
              ("Rh",45);
              ("Pd",46);
              ("Ag",47);
              ("Cd",48);
              ("In",49);
              ("Sn",50);
              ("Sb",51);
              ("Te",52);
              ("I",53);
              ("Xe",54);
              ("Cs",55);
              ("Ba",56);
              ("La",57);
              ("Ce",58);
              ("Pr",59);
              ("Nd",60);
              ("Pm",61);
              ("Sm",62);
              ("Eu",63);
              ("Gd",64);
              ("Tb",65);
              ("Dy",66);
              ("Ho",67);
              ("Er",68);
              ("Tm",69);
              ("Yb",70);
              ("Lu",71);
              ("Hf",72);
              ("Ta",73);
              ("W",74);
              ("Re",75);
              ("Os",76);
              ("Ir",77);
              ("Pt",78);
              ("Au",79);
              ("Hg",80);
              ("Tl",81);
              ("Pb",82);
              ("Bi",83);
              ("Po",84);
              ("At",85);
              ("Rn",86);
              ("Fr",87);
              ("Ra",88);
              ("Ac",89);
              ("Th",90);
              ("Pa",91);
              ("U",92);
              ("Np",93);
              ("Pu",94);
              ("Am",95);
              ("Cm",96);
              ("Bk",97);
              ("Cf",98);
              ("Es",99);
              ("Fm",100);
              ("Md",101);
              ("No",102);
              ("Lr",103);
              ("Rf",104);
              ("Db",105);
              ("Sg",106);
              ("Bh",107);
              ("Hs",108);
              ("Mt",109);
              ("Ds",110);
              ("Rg",111);
              ("Cn",112);
              ("Nh",113);
              ("Fl",114);
              ("Mc",115);
              ("Lv",116);
              ("Ts",117);
              ("Og",118)]

let anum_of_symbol s =
  try Ht.find symbol2anum s
  with Not_found ->
    (Log.fatal "Ptable.anum_of_symbol: no such chemical element: %s" s;
     exit 1)
