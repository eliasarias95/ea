"""
common imports 
"""
from os import path, access, R_OK  # W_OK for write permission.
from operator import itemgetter
import sys

from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from javax.swing import *
from java.awt.Color import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.awt.ColorMap import *

from dnp import LocalSlopeFinder
from utils import Synthetic

