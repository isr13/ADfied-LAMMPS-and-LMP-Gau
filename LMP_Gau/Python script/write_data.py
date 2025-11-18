import os
import convertkit
import sys

inputxyz = sys.argv[1]

mol = convertkit.convertkit(inputxyz)

mol.writeLammpsDataFile()
