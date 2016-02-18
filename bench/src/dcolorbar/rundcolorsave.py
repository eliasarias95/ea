#############################################################################
#Builds PP and PS figures for CWP presentation.


from imports import *
from dcolorbar import Utility
from dcolorbar import SOM1 
from dcolorbar import SOM2

path = "/Users/Chris/Home/box/cg/bench/src/dcolorbar/data/"
def main(args):
  dcolor()

def dcolor():
  #read in color bars
  colorBarNames = ["prism.bin","gray_yellow_red.bin","red_white_blue.bin"]
  nColorBar = len(colorBarNames)
  colorBarRGBValues = zerofloat(3,256,nColorBar)
  colorBarPacked = zerofloat(256*3,nColorBar)
  sumDistanceCB = zerofloat(nColorBar)
  for icb in range(0,nColorBar):
    colorBarRGBValues[icb] = Utility.readL(3,256,path+colorBarNames[icb])
    colorBarPacked[icb] = Utility.readL(256*3,path+colorBarNames[icb])
    rCB,gCB,bCB = getRGB(colorBarRGBValues[icb])
    #test(r,g,b)

    #read in image
    imageFile = "red_white_blue_crop.png"
    #imageFile = "gray_yellow_red_crop.png"
    #imageFile = "prism.png"
    imageRGBValues = Utility.rgbFromFile(path+imageFile)
    rI = Utility.make1D(imageRGBValues[0])
    gI = Utility.make1D(imageRGBValues[1])
    bI = Utility.make1D(imageRGBValues[2])
    rIRand,gIRand,bIRand = randomImage(rI,gI,bI,1000)
    nI = len(rIRand)
    #nI = len(rI)
    nrgbI = 3*nI
    rgbI = zerofloat(nrgbI)
    j = 0
    for i in range(0,nI):
      rgbI[j] = rIRand[i]
      rgbI[j+1] = gIRand[i]
      rgbI[j+2] = bIRand[i]
      #rgbI[j] = rI[i]
      #rgbI[j+1] = gI[i]
      #rgbI[j+2] = bI[i]
      j += 3

    sumDistanceCB[icb] = findSumDistance(rCB,gCB,bCB,rIRand,gIRand,bIRand)
  minInx = minIndex(sumDistanceCB)
  Utility.plot(rgbI,rgbI)
  Utility.plot(colorBarPacked[0],colorBarPacked[0])
  Utility.plot(colorBarPacked[1],colorBarPacked[1])
  Utility.plot(colorBarPacked[2],colorBarPacked[2])
  Utility.plot(colorBarPacked[minInx],colorBarPacked[minInx])
  """
  #send random points into SOM
  niter = 500000
  nh = 8 
  nv = 8
  nattr = 3
  
  som = SOM2(niter,nh,nv,nattr)
  som.train2D(imageRGBValues)
  nodes = som.getNodes()

  weightr = zerofloat(nv,nh)
  weightg = zerofloat(nv,nh)
  weightb = zerofloat(nv,nh)
  for i2 in range(0,nv):
    for i1 in range(0,nh):
      weightr[i1][i2] = nodes[i2][i1][0]
      weightg[i1][i2] = nodes[i2][i1][1]
      weightb[i1][i2] = nodes[i2][i1][2]

  nPacked = 3*nh*nv
  weightRGB = zerofloat(nPacked)
  j = 0
  for ih in range(0,nh):
    for iv in range(0,nv):
      weightRGB[j  ] = weightr[iv][ih]
      weightRGB[j+1] = weightr[iv][ih]
      weightRGB[j+2] = weightr[iv][ih]
      j += 3




  weightr = zerofloat(nh)
  weightg = zerofloat(nh)
  weightb = zerofloat(nh)
  for i in range(0,nh):
    weightr[i] = nodes[i][0]
    weightg[i] = nodes[i][1]
    weightb[i] = nodes[i][2]
    """
def findDistance(r1,g1,b1,r2,g2,b2):
  return sqrt((r1-r2)*(r1-r2)+(g1-g2)*(g1-g2)+(b1-b2)*(b1-b2))

def findSumDistance(rCB,gCB,bCB,rI,gI,bI):
  nCB = len(rCB)
  nI = len(rI)
  sumDistance = 0
  distanceAr = zerofloat(nI,nCB)
  for iI in range(0,nI):
    for iCB in range(0,nCB):
      distanceAr[iCB][iI] = findDistance(rCB[iCB],gCB[iCB],bCB[iCB],rI[iI],gI[iI],bI[iI])
    sumDistance += float(min(distanceAr))
  return sumDistance

def minIndex(array):
  narray = len(array)
  min = FLT_MAX
  index = 0
  for i in range(0,narray):
    if (array[i] < min):
      min = array[i]
      index = i
  return index

def randomImage(rI,gI,bI,nRand):
  nrI = len(rI)
  randRImage = zerofloat(nRand)
  randGImage = zerofloat(nRand)
  randBImage = zerofloat(nRand)
  r = Random()
  for i in range(0,nRand):
    index = r.nextInt(nrI)
    print "index = "+str(index)
    randRImage[i] = rI[index]
    randGImage[i] = gI[index] 
    randBImage[i] = bI[index]
  return randRImage,randGImage,randBImage

def getRGB(colorBarRGBValues):
  nrbg = 3
  nvalues = len(colorBarRGBValues)
  r = zerofloat(nvalues)
  g = zerofloat(nvalues)
  b = zerofloat(nvalues)
  for i in range(0,nvalues):
    r[i] = colorBarRGBValues[i][0]
    g[i] = colorBarRGBValues[i][1]
    b[i] = colorBarRGBValues[i][2]
  return r,g,b

def test(r,g,b):
  cm = ColorMap(0.0,1.0,r,g,b)
  icm = cm.getColorModel()
  rampFloats0To1 = rampfloat(0.0,0.01,100)
  ramp2D = zerofloat(100,2)
  for i in range(0,2):
    ramp2D[i] = rampFloats0To1
  
  pv = PixelsView(ramp2D)
  pv.setColorModel(icm)
  pp = PlotPanel()
  pp.addTiledView(pv)
  pf = PlotFrame(pp)
  pf.setVisible(True)



  



#############################################################################
# Do everything on Swing thread.


class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
