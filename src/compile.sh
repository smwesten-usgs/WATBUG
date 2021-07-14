#!/bin/bash
gfortran watbug__modified.f -fdefault-real-8 -ffixed-form -ffixed-line-length-none -std=legacy -fmax-errors=4 -o watbug
