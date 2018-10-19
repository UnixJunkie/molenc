(* SYBYL atom types
   from the MOL2 file format specification p 53-54 SYBYL 7.1 (Mid-2005)
   https://github.com/UnixJunkie/mol2-file-format-spec *)

module SMap = BatMap.String

type t =
  | C3 (* sp3 carbon *)
  | C2 (* sp2 carbon *)
  | C1 (* sp carbon *)
  | Car (* aromatic carbon *)
  | Ccat (* carbocation used only in a guadinium group *)
  | N3 (* sp3 nitrogen *)
  | N2 (* sp2 nitrogen *)
  | N1 (* sp nitrogen *)
  | Nar (* aromatic nitrogen *)
  | Nam (* amide nitrogen *)
  | Npl3 (* trigonal planar nitrogen *)
  | N4 (* sp3 nitrogen positively charged *)
  | O3 (* sp3 oxygen *)
  | O2 (* sp2 oxygen *)
  | Oco2 (* oxygen in carboxylate and phosphate groups *)
  | Ospc (* oxygen in single point charge water model *)
  | Ot3p (* oxygen in transferable intermolecular potential
            (TIP3P) water model *)
  | S3 (* sp3 sulfur *)
  | S2 (* sp2 sulfur *)
  | SO (* sulfoxide sulfur *)
  | SO2 (* sulfone sulfur *)
  | P3 (* phosphorous sp3 *)
  | F (* fluorine *)
  | H (* hydrogen *)
  | Hspc
  | Ht3p
  | LP (* lone pair *)
  | Du (* dummy atom *)
  | DuC (* dummy carbon *)
  | Any (* any atom *)
  | Hal (* halogen *)
  | Het (* heteroatom: N or O or S or P *)
  | Li (* lithium *)
  | Na (* Sodium *)
  | Mg (* magnesium *)
  | Al (* aluminium *)
  | Si (* silicon *)
  | K (* potassium *)
  | Ca (* calcium *)
  | Crth (* tetrahedral chromium *)
  | Croh (* octahedral chromium *)
  | Mn (* manganese *)
  | Fe (* iron *)
  | Cooh (* octahedral cobalt *)
  | Cu (* copper *)
  | Cl (* chlorine *)
  | Br (* bromine *)
  | I (* iodine *)
  | Zn (* zinc *)
  | Se (* selenium *)
  | Mo (* molybdenum *)
  | Sn (* tin *)
  (* the following atom types are not in the MOL2 spec *)
  | Hg (* mercury *)
  | Nd (* Neodymium *)
  | As (* Arsenic *)
  | B (* Boron *)
  | Ba (* Barium *)
  | Pt (* Platinum *)
  | Ni (* Nickel *)
  | Co (* Cobalt *)
  | Pd (* Palladium *)
  | Au (* Gold *)
  | Ag (* Silver *)
  | Bi (* Bismuth *)
  | Ce (* Cerium *)
  | Cr (* Chrome *)
  | Cs (* Cesium *)
  | Dy (* Dysprosium *)
  | Ga (* Gallium *)
  | Gd (* Gadolinium *)
  | Ge (* Germanium *)
  | He (* Helium *)
  | In (* Indium *)
  | Ir (* Iridium *)
  | Kr (* Krypton *)
  | La (* Lanthanum *)
  | Lu (* Lutetium *)
  | Pb (* Lead *)
  | Ra (* Radium *)
  | Rb (* Rubidium *)
  | Sb (* Antimony *)
  | Sm (* Samarium *)
  | Sr (* Strontium *)
  | Tc (* Technetium *)
  | Ti (* Titanium *)
  | Tl (* Thallium *)
  | Xe (* Xenon *)
  | Zr (* Zirconium *)
  | S (* Sulfur *)
  | V (* Vanadium *)
  | Ac (* Actinium *)
  | Am (* Americium *)
  | Ar (* Argon *)
  | At (* Astatine *)
  | Be (* Beryllium *)
  | Bh (* Bohrium *)
  | Bk (* Berkelium *)
  | Cd (* Cadmium *)
  | Cf (* Californium *)
  | Cm (* Curium *)
  | Db (* Dubnium *)
  | Er (* Erbium *)
  | Es (* Einsteinium *)
  | Eu (* Europium *)
  | Fm (* Fermium *)
  | Fr (* Francium *)
  | Hf (* Hafnium *)
  | Ho (* Holmium *)
  | Hs (* Hassium *)
  | Lr (* Lawrencium *)
  | Md (* Mendelevium *)
  | Mt (* Meitnerium *)
  | Nb (* Niobium *)
  | Ne (* Neon *)
  | No (* Nobelium *)
  | Np (* Neptunium *)
  | Os (* Osmium *)
  | Pa (* Protactinium *)
  | Pm (* Promethium *)
  | Po (* Polonium *)
  | Pr (* Praseodymium *)
  | Pu (* Plutonium *)
  | Re (* Rhenium *)
  | Rf (* Rutherfordium *)
  | Rh (* Rhodium *)
  | Rn (* Radon *)
  | Ru (* Ruthenium *)
  | Sc (* Scandium *)
  | Sg (* Seaborgium *)
  | Ta (* Tantalum *)
  | Tb (* Terbium *)
  | Te (* Tellurium *)
  | Th (* Thorium *)
  | Tm (* Thulium *)
  | U (* Uranium *)
  | W (* Tungsten *)
  | Y (* Yttrium *)
  | Yb (* Ytterbium *)

exception Unknown_atom of string

let of_string = function
  | "C.3"   -> C3
  | "C.2"   -> C2
  | "C.1"   -> C1
  | "C.ar"  -> Car
  | "C.cat" -> Ccat
  | "N.3"   -> N3
  | "N.2"   -> N2
  | "N.1"   -> N1
  | "N.ar"  -> Nar
  | "N.am"  -> Nam
  | "N.pl3" -> Npl3
  | "N.4"   -> N4
  | "O.3"   -> O3
  | "O.2"   -> O2
  | "O.co2" -> Oco2
  | "O.spc" -> Ospc
  | "O.t3p" -> Ot3p
  | "S.3"   -> S3
  | "S.2"   -> S2
  | "S.o"
  | "S.O"   -> SO
  | "S.o2"
  | "S.O2"  -> SO2
  | "P.3"   -> P3
  | "Pt"    -> Pt
  | "F"     -> F
  | "H"     -> H
  | "Hg"    -> Hg
  | "H.spc" -> Hspc
  | "H.t3p" -> Ht3p
  | "LP"    -> LP
  | "Du"    -> Du
  | "Du.C"  -> DuC
  | "Any"   -> Any
  | "Hal"   -> Hal
  | "Het"   -> Het
  | "Li"    -> Li
  | "Na"    -> Na
  | "Nd"    -> Nd
  | "Mg"    -> Mg
  | "Al"    -> Al
  | "As"    -> As
  | "Si"    -> Si
  | "K"     -> K
  | "Ca"    -> Ca
  | "Cr.th" -> Crth
  | "Cr.oh" -> Croh
  | "Mn"    -> Mn
  | "Fe"    -> Fe
  | "Co.oh" -> Cooh
  | "Cu"    -> Cu
  | "Cl"    -> Cl
  | "B"     -> B
  | "Br"    -> Br
  | "Ba"    -> Ba
  | "I"     -> I
  | "Zn"    -> Zn
  | "Se"    -> Se
  | "Mo"    -> Mo
  | "Sn"    -> Sn
  | "Ni"    -> Ni
  | "Co"    -> Co
  | "Pd"    -> Pd
  | "Au"    -> Au
  | "Ag"    -> Ag
  | "Bi"    -> Bi
  | "Ce"    -> Ce
  | "Cr"    -> Cr
  | "Cs"    -> Cs
  | "Dy"    -> Dy
  | "Ga"    -> Ga
  | "Gd"    -> Gd
  | "Ge"    -> Ge
  | "He"    -> He
  | "In"    -> In
  | "Ir"    -> Ir
  | "Kr"    -> Kr
  | "La"    -> La
  | "Lu"    -> Lu
  | "Pb"    -> Pb
  | "Ra"    -> Ra
  | "Rb"    -> Rb
  | "Sb"    -> Sb
  | "Sm"    -> Sm
  | "Sr"    -> Sr
  | "Tc"    -> Tc
  | "Ti"    -> Ti
  | "Tl"    -> Tl
  | "Xe"    -> Xe
  | "Zr"    -> Zr
  | "S"     -> S
  | "V"     -> V
  | "Ac"    -> Ac
  | "Am"    -> Am
  | "Ar"    -> Ar
  | "At"    -> At
  | "Be"    -> Be
  | "Bh"    -> Bh
  | "Bk"    -> Bk
  | "Cd"    -> Cd
  | "Cf"    -> Cf
  | "Cm"    -> Cm
  | "Db"    -> Db
  | "Er"    -> Er
  | "Es"    -> Es
  | "Eu"    -> Eu
  | "Fm"    -> Fm
  | "Fr"    -> Fr
  | "Hf"    -> Hf
  | "Ho"    -> Ho
  | "Hs"    -> Hs
  | "Lr"    -> Lr
  | "Md"    -> Md
  | "Mt"    -> Mt
  | "Nb"    -> Nb
  | "Ne"    -> Ne
  | "No"    -> No
  | "Np"    -> Np
  | "Os"    -> Os
  | "Pa"    -> Pa
  | "Pm"    -> Pm
  | "Po"    -> Po
  | "Pr"    -> Pr
  | "Pu"    -> Pu
  | "Re"    -> Re
  | "Rf"    -> Rf
  | "Rh"    -> Rh
  | "Rn"    -> Rn
  | "Ru"    -> Ru
  | "Sc"    -> Sc
  | "Sg"    -> Sg
  | "Ta"    -> Ta
  | "Tb"    -> Tb
  | "Te"    -> Te
  | "Th"    -> Th
  | "Tm"    -> Tm
  | "U"     -> U
  | "W"     -> W
  | "Y"     -> Y
  | "Yb"    -> Yb
  | unk -> raise (Unknown_atom unk)
  (* (\* code in case you need to catch many at once *\)
   * | _ -> (Log.error "Sybyl: unknown atom: %s" unk; Any) *)

let to_string = function
  | C3   -> "C.3"
  | C2   -> "C.2"
  | C1   -> "C.1"
  | Car  -> "C.ar"
  | Ccat -> "C.cat"
  | N3   -> "N.3"
  | N2   -> "N.2"
  | N1   -> "N.1"
  | Nar  -> "N.ar"
  | Nam  -> "N.am"
  | Npl3 -> "N.pl3"
  | N4   -> "N.4"
  | O3   -> "O.3"
  | O2   -> "O.2"
  | Oco2 -> "O.co2"
  | Ospc -> "O.spc"
  | Ot3p -> "O.t3p"
  | S3   -> "S.3"
  | S2   -> "S.2"
  | SO   -> "S.O"
  | SO2  -> "S.O2"
  | P3   -> "P.3"
  | Pt   -> "Pt"
  | F    -> "F"
  | H    -> "H"
  | Hg   -> "Hg"
  | Hspc -> "H.spc"
  | Ht3p -> "H.t3p"
  | LP   -> "LP"
  | Du   -> "Du"
  | DuC  -> "Du.C"
  | Any  -> "Any"
  | Hal  -> "Hal"
  | Het  -> "Het"
  | Li   -> "Li"
  | Na   -> "Na"
  | Nd   -> "Nd"
  | Mg   -> "Mg"
  | Al   -> "Al"
  | As   -> "As"
  | Si   -> "Si"
  | K    -> "K"
  | Ca   -> "Ca"
  | Crth -> "Cr.th"
  | Croh -> "Cr.oh"
  | Mn   -> "Mn"
  | Fe   -> "Fe"
  | Cooh -> "Co.oh"
  | Cu   -> "Cu"
  | Cl   -> "Cl"
  | B    -> "B"
  | Ba   -> "Ba"
  | Br   -> "Br"
  | I    -> "I"
  | Zn   -> "Zn"
  | Se   -> "Se"
  | Mo   -> "Mo"
  | Sn   -> "Sn"
  | Ni   -> "Ni"
  | Co   -> "Co"
  | Pd   -> "Pd"
  | Au   -> "Au"
  | Ag   -> "Ag"
  | Bi    -> "Bi"
  | Ce    -> "Ce"
  | Cr    -> "Cr"
  | Cs    -> "Cs"
  | Dy    -> "Dy"
  | Ga    -> "Ga"
  | Gd    -> "Gd"
  | Ge    -> "Ge"
  | He    -> "He"
  | In    -> "In"
  | Ir    -> "Ir"
  | Kr    -> "Kr"
  | La    -> "La"
  | Lu    -> "Lu"
  | Pb    -> "Pb"
  | Ra    -> "Ra"
  | Rb    -> "Rb"
  | Sb    -> "Sb"
  | Sm    -> "Sm"
  | Sr    -> "Sr"
  | Tc    -> "Tc"
  | Ti    -> "Ti"
  | Tl    -> "Tl"
  | Xe    -> "Xe"
  | Zr    -> "Zr"
  | S     -> "S"
  | V     -> "V"
  | Ac    -> "Ac"
  | Am    -> "Am"
  | Ar    -> "Ar"
  | At    -> "At"
  | Be    -> "Be"
  | Bh    -> "Bh"
  | Bk    -> "Bk"
  | Cd    -> "Cd"
  | Cf    -> "Cf"
  | Cm    -> "Cm"
  | Db    -> "Db"
  | Er    -> "Er"
  | Es    -> "Es"
  | Eu    -> "Eu"
  | Fm    -> "Fm"
  | Fr    -> "Fr"
  | Hf    -> "Hf"
  | Ho    -> "Ho"
  | Hs    -> "Hs"
  | Lr    -> "Lr"
  | Md    -> "Md"
  | Mt    -> "Mt"
  | Nb    -> "Nb"
  | Ne    -> "Ne"
  | No    -> "No"
  | Np    -> "Np"
  | Os    -> "Os"
  | Pa    -> "Pa"
  | Pm    -> "Pm"
  | Po    -> "Po"
  | Pr    -> "Pr"
  | Pu    -> "Pu"
  | Re    -> "Re"
  | Rf    -> "Rf"
  | Rh    -> "Rh"
  | Rn    -> "Rn"
  | Ru    -> "Ru"
  | Sc    -> "Sc"
  | Sg    -> "Sg"
  | Ta    -> "Ta"
  | Tb    -> "Tb"
  | Te    -> "Te"
  | Th    -> "Th"
  | Tm    -> "Tm"
  | U     -> "U"
  | W     -> "W"
  | Y     -> "Y"
  | Yb    -> "Yb"

let to_int = function
  | C3   -> 0
  | C2   -> 1
  | C1   -> 2
  | Car  -> 3
  | Ccat -> 4
  | N3   -> 5
  | N2   -> 6
  | N1   -> 7
  | Nar  -> 8
  | Nam  -> 9
  | Npl3 -> 10
  | N4   -> 11
  | O3   -> 12
  | O2   -> 13
  | Oco2 -> 14
  | Ospc -> 15
  | Ot3p -> 16
  | S3   -> 17
  | S2   -> 18
  | SO   -> 19
  | SO2  -> 20
  | P3   -> 21
  | Pt   -> 22
  | F    -> 23
  | H    -> 24
  | Hg   -> 25
  | Hspc -> 26
  | Ht3p -> 27
  | LP   -> 28
  | Du   -> 29
  | DuC  -> 30
  | Any  -> 31
  | Hal  -> 32
  | Het  -> 33
  | Li   -> 34
  | Na   -> 35
  | Nd   -> 36
  | Mg   -> 37
  | Al   -> 38
  | As   -> 39
  | Si   -> 40
  | K    -> 41
  | Ca   -> 42
  | Crth -> 43
  | Croh -> 44
  | Mn   -> 45
  | Fe   -> 46
  | Cooh -> 47
  | Cu   -> 48
  | Cl   -> 49
  | B    -> 50
  | Ba   -> 51
  | Br   -> 52
  | I    -> 53
  | Zn   -> 54
  | Se   -> 55
  | Mo   -> 56
  | Sn   -> 57
  | Ni   -> 58
  | Co   -> 59
  | Pd   -> 60
  | Au   -> 61
  | Ag   -> 62
  | Bi   -> 63
  | Ce   -> 64
  | Cr   -> 65
  | Cs   -> 66
  | Dy   -> 67
  | Ga   -> 68
  | Gd   -> 69
  | Ge   -> 70
  | He   -> 71
  | In   -> 72
  | Ir   -> 73
  | Kr   -> 74
  | La   -> 75
  | Lu   -> 76
  | Pb   -> 77
  | Ra   -> 78
  | Rb   -> 79
  | Sb   -> 80
  | Sm   -> 81
  | Sr   -> 82
  | Tc   -> 83
  | Ti   -> 84
  | Tl   -> 85
  | Xe   -> 86
  | Zr   -> 87
  | S    -> 88
  | V    -> 89
  | Ac   -> 90 
  | Am   -> 91 
  | Ar   -> 92 
  | At   -> 93 
  | Be   -> 94 
  | Bh   -> 95 
  | Bk   -> 96 
  | Cd   -> 97 
  | Cf   -> 98 
  | Cm   -> 99 
  | Db   -> 100
  | Er   -> 101
  | Es   -> 102
  | Eu   -> 103
  | Fm   -> 104
  | Fr   -> 105
  | Hf   -> 106
  | Ho   -> 107
  | Hs   -> 108
  | Lr   -> 109
  | Md   -> 110
  | Mt   -> 111
  | Nb   -> 112
  | Ne   -> 113
  | No   -> 114
  | Np   -> 115
  | Os   -> 116
  | Pa   -> 117
  | Pm   -> 118
  | Po   -> 119
  | Pr   -> 120
  | Pu   -> 121
  | Re   -> 122
  | Rf   -> 123
  | Rh   -> 124
  | Rn   -> 125
  | Ru   -> 126
  | Sc   -> 127
  | Sg   -> 128
  | Ta   -> 129
  | Tb   -> 130
  | Te   -> 131
  | Th   -> 132
  | Tm   -> 133
  | U    -> 134
  | W    -> 135
  | Y    -> 136
  | Yb   -> 137

let base = 138

(*
let symbol_to_atomic_number =
  SMap.add "H" 1 SMap.empty |>
  SMap.add "He" 2 |>
  SMap.add "Li" 3 |>
  SMap.add "Be" 4 |>
  SMap.add "B" 5 |>
  SMap.add "C" 6 |>
  SMap.add "N" 7 |>
  SMap.add "O" 8 |>
  SMap.add "F" 9 |>
  SMap.add "Ne" 10 |>
  SMap.add "Na" 11 |>
  SMap.add "Mg" 12 |>
  SMap.add "Al" 13 |>
  SMap.add "Si" 14 |>
  SMap.add "P" 15 |>
  SMap.add "S" 16 |>
  SMap.add "Cl" 17 |>
  SMap.add "Ar" 18 |>
  SMap.add "K" 19 |>
  SMap.add "Ca" 20 |>
  SMap.add "Sc" 21 |>
  SMap.add "Ti" 22 |>
  SMap.add "V" 23 |>
  SMap.add "Cr" 24 |>
  SMap.add "Mn" 25 |>
  SMap.add "Fe" 26 |>
  SMap.add "Co" 27 |>
  SMap.add "Ni" 28 |>
  SMap.add "Cu" 29 |>
  SMap.add "Zn" 30 |>
  SMap.add "Ga" 31 |>
  SMap.add "Ge" 32 |>
  SMap.add "As" 33 |>
  SMap.add "Se" 34 |>
  SMap.add "Br" 35 |>
  SMap.add "Kr" 36 |>
  SMap.add "Rb" 37 |>
  SMap.add "Sr" 38 |>
  SMap.add "Y" 39 |>
  SMap.add "Zr" 40 |>
  SMap.add "Nb" 41 |>
  SMap.add "Mo" 42 |>
  SMap.add "Tc" 43 |>
  SMap.add "Ru" 44 |>
  SMap.add "Rh" 45 |>
  SMap.add "Pd" 46 |>
  SMap.add "Ag" 47 |>
  SMap.add "Cd" 48 |>
  SMap.add "In" 49 |>
  SMap.add "Sn" 50 |>
  SMap.add "Sb" 51 |>
  SMap.add "Te" 52 |>
  SMap.add "I" 53 |>
  SMap.add "Xe" 54 |>
  SMap.add "Cs" 55 |>
  SMap.add "Ba" 56 |>
  SMap.add "La" 57 |>
  SMap.add "Ce" 58 |>
  SMap.add "Pr" 59 |>
  SMap.add "Nd" 60 |>
  SMap.add "Pm" 61 |>
  SMap.add "Sm" 62 |>
  SMap.add "Eu" 63 |>
  SMap.add "Gd" 64 |>
  SMap.add "Tb" 65 |>
  SMap.add "Dy" 66 |>
  SMap.add "Ho" 67 |>
  SMap.add "Er" 68 |>
  SMap.add "Tm" 69 |>
  SMap.add "Yb" 70 |>
  SMap.add "Lu" 71 |>
  SMap.add "Hf" 72 |>
  SMap.add "Ta" 73 |>
  SMap.add "W" 74 |>
  SMap.add "Re" 75 |>
  SMap.add "Os" 76 |>
  SMap.add "Ir" 77 |>
  SMap.add "Pt" 78 |>
  SMap.add "Au" 79 |>
  SMap.add "Hg" 80 |>
  SMap.add "Tl" 81 |>
  SMap.add "Pb" 82 |>
  SMap.add "Bi" 83 |>
  SMap.add "Po" 84 |>
  SMap.add "At" 85 |>
  SMap.add "Rn" 86 |>
  SMap.add "Fr" 87 |>
  SMap.add "Ra" 88 |>
  SMap.add "Ac" 89 |>
  SMap.add "Th" 90 |>
  SMap.add "Pa" 91 |>
  SMap.add "U" 92 |>
  SMap.add "Np" 93 |>
  SMap.add "Pu" 94 |>
  SMap.add "Am" 95 |>
  SMap.add "Cm" 96 |>
  SMap.add "Bk" 97 |>
  SMap.add "Cf" 98 |>
  SMap.add "Es" 99 |>
  SMap.add "Fm" 100 |>
  SMap.add "Md" 101 |>
  SMap.add "No" 102 |>
  SMap.add "Lr" 103 |>
  SMap.add "Rf" 104 |>
  SMap.add "Db" 105 |>
  SMap.add "Sg" 106 |>
  SMap.add "Bh" 107 |>
  SMap.add "Hs" 108 |>
  SMap.add "Mt" 109

let max_AN = 109

let sym2AN = symbol_to_atomic_number

let atomic_number = function
  | C3 (* sp3 carbon *)
  | C2 (* sp2 carbon *)
  | C1 (* sp carbon *)
  | Car (* aromatic carbon *)
  | Ccat (* carbocation used only in a guadinium group *)
    -> SMap.find "C" sym2AN
  | N3 (* sp3 nitrogen *)
  | N2 (* sp2 nitrogen *)
  | N1 (* sp nitrogen *)
  | Nar (* aromatic nitrogen *)
  | Nam (* amide nitrogen *)
  | Npl3 (* trigonal planar nitrogen *)
  | N4 (* sp3 nitrogen positively charged *)
    -> SMap.find "N" sym2AN
  | O3 (* sp3 oxygen *)
  | O2 (* sp2 oxygen *)
  | Oco2 (* oxygen in carboxylate and phosphate groups *)
  | Ospc (* oxygen in single point charge water model *)
  | Ot3p (* oxygen in transferable intermolecular potential
            (TIP3P) water model *)
    -> SMap.find "O" sym2AN
  | S3 (* sp3 sulfur *)
  | S2 (* sp2 sulfur *)
  | SO (* sulfoxide sulfur *)
  | SO2 (* sulfone sulfur *)
    -> SMap.find "S" sym2AN
  | P3 (* phosphorous sp3 *)
    -> SMap.find "P" sym2AN
  | F (* fluorine *)
    -> SMap.find "F" sym2AN
  | H (* hydrogen *)
  | Hspc
  | Ht3p
    -> SMap.find "H" sym2AN
  | LP (* lone pair *) -> failwith "atomic_number of lone pair"
  | Du (* dummy atom *) -> failwith "atomic number of dummy atom"
  | DuC (* dummy carbon *) -> SMap.find "C" sym2AN
  | Any (* any atom *) -> failwith "atomic number of any atom"
  | Hal (* halogen *) -> failwith "atomic number of halogen"
  | Het (* heteroatom: N or O or S or P *)
    -> failwith "atomic number of hetereoatom"
  | Li (* lithium *) -> SMap.find "Li" sym2AN
  | Na (* Sodium *) -> SMap.find "Na" sym2AN
  | Mg (* magnesium *) -> SMap.find "Mg" sym2AN
  | Al (* aluminium *) -> SMap.find "Al" sym2AN
  | Si (* silicon *) -> SMap.find "Si" sym2AN
  | K (* potassium *) -> SMap.find "K" sym2AN
  | Ca (* calcium *) -> SMap.find "Ca" sym2AN
  | Crth (* tetrahedral chromium *)
  | Croh (* octahedral chromium *)
    -> SMap.find "Cr" sym2AN
  | Mn (* manganese *) -> SMap.find "Mn" sym2AN
  | Fe (* iron *) -> SMap.find "Fe" sym2AN
  | Cooh (* octahedral cobalt *) -> SMap.find "Co" sym2AN
  | Cu (* copper *) -> SMap.find "Cu" sym2AN
  | Cl (* chlorine *) -> SMap.find "Cl" sym2AN
  | Br (* bromine *) -> SMap.find "Br" sym2AN
  | I (* iodine *) -> SMap.find "I" sym2AN
  | Zn (* zinc *) -> SMap.find "Zn" sym2AN
  | Se (* selenium *) -> SMap.find "Se" sym2AN
  | Mo (* molybdenum *) -> SMap.find "Mo" sym2AN
  | Sn (* tin *) -> SMap.find "Sn" sym2AN
  (* the following atom types are not in the MOL2 spec *)
  | Hg (* mercury *) -> SMap.find "Hg" sym2AN
  | Nd (* Neodymium *) -> SMap.find "Nd" sym2AN
  | As (* Arsenic *) -> SMap.find "As" sym2AN
  | B (* Boron *) -> SMap.find "B" sym2AN
  | Ba (* Barium *) -> SMap.find "Ba" sym2AN
  | Pt (* Platinum *) -> SMap.find "Pt" sym2AN
  | Ni (* Nickel *) -> SMap.find "Ni" sym2AN
  | Co (* Cobalt *) -> SMap.find "Co" sym2AN
  | Pd (* Palladium *) -> SMap.find "Pd" sym2AN
  | Au (* Gold *) -> SMap.find "Au" sym2AN
  | Ag (* Silver *) -> SMap.find "Ag" sym2AN
  | Bi (* Bismuth *) -> SMap.find "Bi" sym2AN
  | Ce (* Cerium *) -> SMap.find "Ce" sym2AN
  | Cr (* Chrome *) -> SMap.find "Cr" sym2AN
  | Cs (* Cesium *) -> SMap.find "Cs" sym2AN
  | Dy (* Dysprosium *) -> SMap.find "Dy" sym2AN
  | Ga (* Gallium *) -> SMap.find "Ga" sym2AN
  | Gd (* Gadolinium *) -> SMap.find "Gd" sym2AN
  | Ge (* Germanium *) -> SMap.find "Ge" sym2AN
  | He (* Helium *) -> SMap.find "He" sym2AN
  | In (* Indium *) -> SMap.find "In" sym2AN
  | Ir (* Iridium *) -> SMap.find "Ir" sym2AN
  | Kr (* Krypton *) -> SMap.find "Kr" sym2AN
  | La (* Lanthanum *) -> SMap.find "La" sym2AN
  | Lu (* Lutetium *) -> SMap.find "Lu" sym2AN
  | Pb (* Lead *) -> SMap.find "Pb" sym2AN
  | Ra (* Radium *) -> SMap.find "Ra" sym2AN
  | Rb (* Rubidium *) -> SMap.find "Rb" sym2AN
  | Sb (* Antimony *) -> SMap.find "Sb" sym2AN
  | Sm (* Samarium *) -> SMap.find "Sm" sym2AN
  | Sr (* Strontium *) -> SMap.find "Sr" sym2AN
  | Tc (* Technetium *) -> SMap.find "Tc" sym2AN
  | Ti (* Titanium *) -> SMap.find "Ti" sym2AN
  | Tl (* Thallium *) -> SMap.find "Tl" sym2AN
  | Xe (* Xenon *) -> SMap.find "Xe" sym2AN
  | Zr (* Zirconium *) -> SMap.find "Zr" sym2AN
  | S (* Sulfur *) -> SMap.find "S" sym2AN
  | V (* Vanadium *) -> SMap.find "V" sym2AN
*)
