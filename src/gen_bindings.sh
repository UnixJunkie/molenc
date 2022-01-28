#!/bin/bash

pyml_bindgen rdkit_wrapper_specs.txt rdkit_wrapper Rdkit_wrapper --caml-module Rdkit > rdkit.ml
