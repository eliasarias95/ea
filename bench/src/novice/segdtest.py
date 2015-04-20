"""
Reads segd files from CSM field camp, assuming these are in Sercel's 
SEG-D Rev 1 format, and writes a dat file containing a 3D array of 
floats. The byte order for floats in the dat file is BIG_ENDIAN. 
Author: Dave Hale, Colorado School of Mines
Version: 2013.10.15
"""
from imports import *

s1 = Sampling(4001,0.001,0.000) # time sampling
s2 = Sampling(540,0.010,0.0) # station sampling
#s3 = Sampling(10,1.0,1160) # first shotpoint is 1000
#s2 = Sampling(277,0.015,-0.030) # channel sampling
#s3 = Sampling(1,1.0,1001.0) # shotpoint station sampling A
segdDir = "/Users/earias/Home/git/ea/bench/src/novice/data/"

#############################################################################
def main(args):
  readAndPlotSegd()
  #readAndPlotDsu()

def readAndPlotDsu():
  segdList = File(segdDir).listFiles() # list of segd files
  nshot = len(segdList) # ignore first 3 shots
  #print segdList
  for segdFile in segdList[10:11]:
    print segdFile
    sl,sp,rpf,rpl,fs = readSegdDsu(segdFile)
    sl,sp = int(sl),int(sp)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    s2 = Sampling(len(fs[0]),0.010,(rpf-sp)*0.010)
    for f in fs:
      gain2(f)
      #g = ppf2(f)
      #gb = ppf2(f,True)
      #gain2(g)
      #gain2(gb)
    plotDsu(s1,s2,fs,title="Shot "+str(sp)+": raw")
  #plotDsuDots(fs)

def readAndPlotSegd():
  segdList = File(segdDir).listFiles() # list of segd files
  #nshot = len(segdList) # ignore first 3 shots
  #print segdList
  for segdFile in segdList:
    print segdFile
    sl,sp,rpf,rpl,f = readSegd(segdFile)
    sl,sp = int(sl),int(sp)
    print "sl =",sl," sp =",sp," rpf =",rpf," rpl =",rpl
    s2 = Sampling(len(f),0.010,(rpf-sp)*0.010)
    #gain2(f)
#    g = ppf2(f)
#    gb = ppf2(f,True)
    #gain2(g)
    #gain2(gb)
    #plot(s1,s2,f,title="Shot "+str(sp)+": raw")
#    plot(s1,s2,g,title="Shot "+str(sp)+": ppf")
#    plot(s1,s2,gb,title="Shot "+str(sp)+": ppfb")
#    plotAmp(s1,s2,f,title="Shot "+str(sp)+": raw")
#    plotAmp(s1,s2,g,title="Shot "+str(sp)+": ppf")
#    plotAmp(s1,s2,gb,title="Shot "+str(sp)+": ppfb")

def ppf2(f,b=False):
  n1,n2 = len(f[0]),len(f)
  m1,m2 = 40,40
  f1=[0.000, 0.120, 0.120, 0.000]
  f2=[0.000,-0.350, 0.350, 0.000]
  g = copy(f)
  if not b:
    ppf = PolygonPassFilter(n1,n2,m1,m2,f1,f2)
  else:
    ppf = PolygonPassFilterB(n1,n2,f1,f2)
  ppf.apply(f,g)
  return g

def tpow2(f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0,0.002,0.0,n1,n2) # time
  t = pow(t,2.0)
  return mul(t,f)

def tpow3(f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  t = rampfloat(s1.first,s1.delta,0.0,n1,n2) # time
  mul(t,t,t) # time squared
  for f3 in f:
    mul(t,f3,f3)

def gain2(f):
  ref = RecursiveExponentialFilter(4000.0)
  for f2 in f: # for all traces, ...
    if max(abs(f2))>0.0:
      g = mul(f2,f2) # square the trace amplitudes
      ref.apply1(g,g) # smooth
      div(f2,sqrt(g),f2) # normalizes

def gain3(f):
  ref = RecursiveExponentialFilter(40.0)
  for f3 in f:
    if max(abs(f3))>0.0:
      g = mul(f3,f3)
      ref.apply1(g,g)
      div(f3,sqrt(g),f3)

def lowpass2(f):
  f3db = 25.0*0.002
  #f3db = 35.0*0.002
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def lowpass3(f):
  bf = ButterworthFilter(0.05,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def plot(s1,s2,f,title=None):
  print "plot f: min =",min(f),"max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  #sp.setSize(750,1000)
  sp.setSize(350,700)
  sp.setVLabel("Time (s)")
  if s2.delta==1.0:
    sp.setHLabel("Station")
  else:
    sp.setHLabel("Offset (km)")
  sp.setVLimits(0.0,1.0)
  #sp.setHLimits(s2.first,s2.first+1.5)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(s1,s2,f)
  #pv.setColorModel(ColorMap.BLUE_WHITE_RED)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setPercentiles(02,98)
  #pv.setClips(-2.5,2.5)

def plotDsuDots(f):
  x1 = rampfloat(1242.0,1.0,100)
  x2 = zerofloat(100)
  sv = Seismic3cView(x1,x2,mul(0.5,f))
  sv.setSample(300)
  pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  pp.setTitle("Particle motion")
  #pp.setVLabel("Northing")
  pp.setHLabel("Station")
  pp.setHLimits(1242,1292)
  pp.setVLimits(-50,50)
  pp.addTiledView(sv)
  pf = PlotFrame(pp)
  pf.setSize(1300,800)
  pf.setVisible(True)
  ss = JSlider(0,1000,1)
  ss.setOrientation(JSlider.HORIZONTAL)
  pf.add(ss,BorderLayout.SOUTH)
  class Changer(ChangeListener):
    def stateChanged(self,e):
      sv.setSample(ss.value)
  ss.addChangeListener(Changer())
  ns = len(f[0][0])
  for ks in range(0,ns):
    sv.setSample(ks)
  for ks in range(ns-1,-1,-1):
    sv.setSample(ks)

def plotDsu(s1,s2,f,title=None):
  print "plot f: min =",min(f),"max =",max(f)
  pp = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  if title:
    pp.setTitle(title+" (vertical, inline, crossline)")
  else:
    pp.setTitle("vertical, inline, crossline")
  pp.setVLabel("Time (s)")
  pp.setHLabel(0,"Offset (km)")
  pp.setHLabel(1,"Offset (km)")
  pp.setHLabel(2,"Offset (km)")
  pp.setVLimits(0.0,1.5)
  pv0 = pp.addPixels(0,0,s1,s2,f[0])
  pv1 = pp.addPixels(0,1,s1,s2,f[1])
  pv2 = pp.addPixels(0,2,s1,s2,f[2])
  pv0.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv1.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv2.setInterpolation(PixelsView.Interpolation.NEAREST)
  #pv.setPercentiles(02,98)
  pv0.setClips(-2.0,2.0)
  pv1.setClips(-2.0,2.0)
  pv2.setClips(-2.0,2.0)
  pf = PlotFrame(pp)
  pf.setSize(1300,800)
  pf.setVisible(True)
  #sp.setHLimits(s2.first,s2.first+1.5)

def plotAmp(s1,s2,f,title=None):
  #s1 = Sampling(s1.count,1.0,0.0)
  #s2 = Sampling(s2.count,1.0,0.0)
  fft = Fft(s1,s2)
  fft.setCenter(True)
  sf1 = fft.getFrequencySampling1()
  sf2 = fft.getFrequencySampling2()
  a = cabs(fft.applyForward(f))
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(700,700)
  sp.setVLabel("Frequency (Hz)")
  sp.setHLabel("Wavenumber (cycles/km)")
  #sp.setVLimits(0.0,60.0)
  if title:
    sp.setTitle(title)
  pv = sp.addPixels(sf1,sf2,a)
  pv.setColorModel(ColorMap.JET)
  pv.setPercentiles(1,99)
  #pv.setClips(-2.5,2.5)

def readData(s1,s2,s3,fileName,bo=ByteOrder.LITTLE_ENDIAN):
  n1,n2,n3 = s1.count,s2.count,s3.count
  f = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName,bo)
  ais.readFloats(f)
  ais.close()
  return f

def writeData(flist,fileName,bo=ByteOrder.LITTLE_ENDIAN):
  n3 = len(flist)
  print "writing",n3," shot records to",fileName
  aos = ArrayOutputStream(fileName,bo)
  for f in flist:
    aos.writeFloats(f)
  aos.close()

def readSegdDsu(segdFile):
  n1,n2 = 4001,300 # number of samples, number of traces
  gh = zerobyte(32) # general header
  th = zerobyte(20) # trace header
  the = zerobyte(32) # trace header extension
  csh = zerobyte(32) # channel set header
  ais = ArrayInputStream(segdFile,ByteOrder.BIG_ENDIAN)
  ais.readBytes(gh) # general header 1
  fn = bcd2(gh,0) # file number
  ais.readBytes(gh) # general header 2
  ais.readBytes(gh) # general header 3
  sln = bin5(gh,3) # source line number
  spn = bin5(gh,8) # source point number
  print "file =",segdFile
  print "fn = ",fn," sln =",sln," spn =",spn
  cns = 2 # channel set number for DSU seismic traces
  ncs = 300 # number of channels = 100 * 3 components
  nct = 542 # total number of channels (including aux channels)
  for ics in range(16): # for each channel set header, ...
    ais.readBytes(csh) # read channel set header
    cn = csh[1] # channel set number
    ct = (csh[10]>>4)&0xf # channel type (in high 4 bits)
    nc = bcd2(csh,8) # number of channels
  ais.skipBytes(1024) # skip extended header
  ais.skipBytes(1024) # skip external header
  rpf = 1
  rpl = 1
  f = None
  for ict in range(nct): # for all channels (including aux channels)
    ais.readBytes(th) # trace header
    cn = th[3] # channel set number
    ic = bcd2(th,4) # channel (trace) number
    ais.readBytes(the) # trace header extension 1
    rln = bin3(the,0) # receiver line number
    rpn = bin3(the,3) # receiver point number
    n1 = bin3(the,7) # number of samples
    #print "ic =",ic," rln =",rln," rpn =",rpn," n1 =",n1
    if ic==1:
      rpf = rpn
    elif ic==n2:
      rpl = rpn
    ais.skipBytes(6*len(the)) # skip trace header extensions 2-7
    if cn==cns: # if seismic DSU channel, ...
      #print "ic =",ic," rln =",rln," rpn =",rpn
      if not f:
        f = zerofloat(n1,n2) # the traces
      ais.readFloats(f[ic-1]) # get the seismic trace
    else:
      ais.skipBytes(4*n1) # skip the auxiliary trace
  ais.close()
  f = mul(1.0e-14,f) # scale values to approximate range [-10,10]
  g = zerofloat(4001,100,3)
  for i in range(100):
    g[0][i] = f[i*3] # vertical
    g[1][i] = f[i*3+1] # ???
    g[2][i] = f[i*3+2] # ???
  return sln,spn,rpf,rpl,g

def readSegd(segdFile):
  n1,n2 = 4001,540 # number of samples, number of traces
  gh = zerobyte(32) # general header
  th = zerobyte(20) # trace header
  the = zerobyte(32) # trace header extension
  csh = zerobyte(32) # channel set header
  ais = ArrayInputStream(segdFile,ByteOrder.BIG_ENDIAN)
  ais.readBytes(gh) # general header 1
  fn = bcd2(gh,0) # file number
  ais.readBytes(gh) # general header 2
  ais.readBytes(gh) # general header 3
  sln = bin5(gh,3) # source line number
  spn = bin5(gh,8) # source point number
  print "file =",segdFile
  print "fn = ",fn," sln =",sln," spn =",spn
  cns = 0 # channel set number for seismic traces
  nct = 0 # total number of channels, including aux channels
  for ics in range(16): # for each channel set header, ...
    ais.readBytes(csh) # read channel set header
    cn = csh[1] # channel set number
    ct = (csh[10]>>4)&0xf # channel type (in high 4 bits)
    nc = bcd2(csh,8) # number of channels
    if nc>0: # if we have channels of this type, ...
      print "cn =",cn," nc =",nc," ct =",ct
      if ct==1: # if seismic, ...
        cns = cn # remember channel set number for seismic
        ncs = nc # remember number of seismic channels
      nct += nc # count total number of channels
  print "nct =",nct,"cns =",cns
  ais.skipBytes(1024) # skip extended header
  ais.skipBytes(1024) # skip external header
  rpf = 1
  rpl = 1
  f = None
  for ict in range(nct): # for all channels (including aux channels)
    ais.readBytes(th) # trace header
    cn = th[3] # channel set number
    ic = bcd2(th,4) # channel (trace) number
    ais.readBytes(the) # trace header extension 1
    rln = bin3(the,0) # receiver line number
    rpn = bin3(the,3) # receiver point number
    n1 = bin3(the,7) # number of samples
    #print "ic =",ic," rln =",rln," rpn =",rpn," n1 =",n1
    if ic==1:
      rpf = rpn
    elif ic==n2:
      rpl = rpn
    ais.skipBytes(6*len(the)) # skip trace header extensions 2-7
    if cn==cns: # if seismic channel, ...
      #print "ic =",ic," rln =",rln," rpn =",rpn
      if not f:
        f = zerofloat(n1,n2) # the traces
      ais.readFloats(f[ic-1]) # get the seismic trace
    else:
      ais.skipBytes(4*n1) # skip the auxiliary trace
  ais.close()
  f = mul(1.0e-14,f) # scale values to approximate range [-10,10]
  return sln,spn,rpf,rpl,f

def bcd2(b,k):
  """ Returns binary-coded decimal integer from bytes k,k+1 in b."""
  return (1000*((b[k  ]>>4)&0xf)+100*(b[k  ]&0xf)+
            10*((b[k+1]>>4)&0xf)+  1*(b[k+1]&0xf))

def displayLine10(vib):
  if vib=="A":
    f = readData(s1,s2,s3a,shotDir+"shota.dat")
  elif vib=="B":
    f = readData(s1,s2,s3b,shotDir+"shotb.dat")
  lowpass3(f)
  tpow3(f)
  gain3(f)
  sf = SimpleFrame()
  ip = sf.addImagePanels(f)
  ip.setPercentiles(2,98)

def bin3(b,k):
  """ Returns integer value encoded in three bytes b[k:k+2]. """
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2]
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  return b0*65536+b1*256+b2

def bin5(b,k):
  """ 
  Returns fixed-point value [0-99999.99] encoded in five bytes b[k:k+4]. 
  """
  b0 = b[k  ]
  b1 = b[k+1]
  b2 = b[k+2] 
  b3 = b[k+3] 
  b4 = b[k+4] 
  if b0<0: b0 += 256
  if b1<0: b1 += 256
  if b2<0: b2 += 256
  if b3<0: b3 += 256
  if b4<0: b4 += 256
  return b0*65536.0+b1*256.0+b2+b3/256.0+b4/65536.0

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
