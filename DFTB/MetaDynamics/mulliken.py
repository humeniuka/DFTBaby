#!/usr/bin/env python

import numpy as np
import os

try:
    from DFTB import utils, DFTB2
except ImportError:
    pass


def get_new_charges_dftb((i, symbols, coords)):
    elements={"H": 1, "C": 6, "N": 7, "O": 8, "Cl": 17, "Br": 35, "Ru": 44}
    atomlist = [(elements[sym.capitalize()], tuple(coord)) for sym, coord in zip(symbols, coords)]
    parser = utils.OptionParserFuncWrapper([DFTB2.DFTB2.__init__], "")
    (options, args) = parser.parse_args()
    dftb = DFTB2.DFTB2(atomlist, **options)
    dftb.setGeometry(atomlist)
    options = {}
    dftb.getEnergy(**options)
    return dftb.dq

