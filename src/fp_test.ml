(* Copyright (C) 2020, Francois Berenger

   Yamanishi laboratory,
   Department of Bioscience and Bioinformatics,
   Faculty of Computer Science and Systems Engineering,
   Kyushu Institute of Technology,
   680-4 Kawazu, Iizuka, Fukuoka, 820-8502, Japan. *)

(* regression tests for Fingerprint.tanimoto *)

module Fp = Molenc.Fingerprint
module Log = Dolog.Log

let () =
  Log.color_on ();
  Log.set_log_level Log.INFO;
  let tani = Fp.tanimoto in
  let fp1 = Fp.of_string "[0:1;1:2;5:1;11:4]" in
  let fp2 = Fp.of_string "[1:1;3:2;11:3]" in
  let fp3 = Fp.of_string "[]" in
  let fp4 = Fp.of_string "[1:2;5:1;11:2]" in
  assert(tani fp1 fp1 = 1.0);
  assert(tani fp2 fp2 = 1.0);
  assert(tani fp3 fp3 = 0.0);
  assert(tani fp4 fp4 = 1.0);
  assert(tani fp1 fp2 = 4.0 /. 10.0);
  assert(tani fp1 fp3 = 0.0);
  assert(tani fp1 fp4 = 5.0 /. 8.0);
  Log.info "OK"
