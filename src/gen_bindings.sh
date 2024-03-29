#!/bin/bash

# pyml_bindgen rdkit_wrapper_specs.txt rdkit_wrapper Rdkit \
#              --caml-module=Rdkit --of-pyo-ret-type=no_check > rdkit.ml

pyml_bindgen rdkit_wrapper_specs.txt rdkit_wrapper Rdkit \
             --embed-python-source rdkit_wrapper.py \
             --caml-module=Rdkit --of-pyo-ret-type=no_check > rdkit.ml

# format the generated code
ocamlformat --inplace --enable-outside-detected-project rdkit.ml
