import os
import convertkit_mm
import sys

inputxyz = sys.argv[1]

mol = convertkit_mm.convertkit(inputxyz)

mol.writeLammpsDataFile()
