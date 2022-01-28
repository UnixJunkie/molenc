#!/bin/bash

pyml_bindgen rdkit_wrapper_specs.txt rdkit_wrapper Rdkit \
             --caml-module=Rdkit --of-pyo-ret-type=no_check > rdkit.ml
