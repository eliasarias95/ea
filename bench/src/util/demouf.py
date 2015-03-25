"""
Real example for unfaulting and unfolding
Author: Xinming Wu, Colorado School of Mines
Version: 2015.02.24
"""
from uff import *
from util import *
from schutils import *
setupForSubset("suf")
s1,s2,s3 = getSamplings()
n1,n2,n3 = s1.count,s2.count,s3.count

# Names and descriptions of files.
g0file  = "g0" # raw input image
gxfile  = "gx" # input image, after bilateral filtering
gsxfile = "gsx" # image after lsf with sharp faults
fwfile  = "fw" # image after unfaulting
fwcfile  = "fwc" # image after unfaulting
fwsfile  = "fws" # image after unfaulting
gwfile  = "gw" # image after unfolding
epfile  = "ep" # eigenvalue-derived planarity
p2file  = "p2" # inline slopes
p3file  = "p3" # crossline slopes
flfile  = "fl" # fault likelihood
fifile  = "fic" # flattened image
gffile  = "gf" # flattened image
wsfile  = "ws" # weight image for flattening
fpfile  = "fp" # fault strike (phi)
ftfile  = "ft" # fault dip (theta)
fltfile = "flt" # fault likelihood thinned
fptfile = "fpt" # fault strike thinned
fttfile = "ftt" # fault dip thinned
fs1file = "fs1" # fault slip (1st component)
fs2file = "fs2" # fault slip (2nd component)
fs3file = "fs3" # fault slip (3rd component)
ft1file = "ft1" # fault slip interpolated (1st component)
ft2file = "ft2" # fault slip interpolated (2nd component)
ft3file = "ft3" # fault slip interpolated (3rd component)
fr1file = "fr1" # fault slip interpolated (1st component)
fr2file = "fr2" # fault slip interpolated (2nd component)
fr3file = "fr3" # fault slip interpolated (3rd component)
fc1file = "fc1"
fc2file = "fc2"
fc3file = "fc3"
fskbase = "fsk" # fault skin (basename only)
fslbase = "fsl" # fault skins after reskinning (basename only)
fskgood = "fsd" # fault skins after reskinning (basename only)
fsksubs = "fss"
fspbase = "fsp"
gxsfile = "gxs"


# These parameters control the scan over fault strikes and dips.
sigmaPhi,sigmaTheta = 20,40
minPhi,maxPhi = 0,360
minTheta,maxTheta = 65,85

# These parameters control the construction of fault skins.
lowerLikelihood = 0.01
upperLikelihood = 0.20
#minSkinSize = 10000
minSkinSize = 200

# These parameters control the computation of fault dip slips.
minThrow =0.0
maxThrow =38.0


# Directory for saved pn images. If None, png images will not be saved;
# otherwise, must create the specified directory before running this script.
#pngDir = None
#pngDir = "../../../png/sch/suf/"
pngDir = "../../../png/sch/suf/slices/"

# We can avoid most computations entirely be setting plotOnly to True.
plotOnly = False

# Processing begins here. When experimenting with one part of this demo, we
# can comment out other parts that have already written results to files.
def main(args):
  #goSubset()
  #goSmooth()
  #goSlip()
  #goUnfaultS()
  #goUnfaultC()
  #goUnfold()
  #goTeaser()
  #goSlices()
  #goShifts()
  goHorizon()
  #goView3D()

def goView3D():
  mark = -100
  gx = readImage(gxfile)
  gw = readImage(gwfile)
  fw = readImage(fwsfile)
  t1 = readImage(ft1file)
  skp = readSkins(fspbase)
  skl = readSkins(fslbase)
  fss = like(gx)
  fss = fillfloat(mark,n1,n2,n3)
  FaultSkin.getThrow(mark,skl,fss)
  d1 = 0.002
  mul(t1,d1*1000,t1)
  plot3(gx,slices=[93,180,192],clab="Amplitude",png="image3D")
  plot3(fw,slices=[93,180,192],png="unfault3D")
  #plot3(gw,slices=[93,190,175],png="unfold3D")
  plot3(gw,slices=[93,180,192],png="unfold3D")
  plot3(gx,t1,cmin=-24,cmax=18,cmap=jetFill(0.3),slices=[93,180,192],
        clab="Vertical shift (ms)",png="shifts3D")
  plot3(gx,skins=skp,smax=20,slices=[93,180,192],png="faultSurfs")
  mul(fss,d1*1000,fss)
  cmin=-0.01
  cmax=max(fss)-20
  plot3(gx,fss,cmin=cmin,cmax=cmax,cmap=jetFillExceptMin(1.0),
        slices=[93,180,192],clab="Fault throw (ms)",cint=5,png="throw3D")

def goHorizon():
  hs=[16,43,63,93,117]
  gx = readImage(gxfile)
  gw = readImage(gwfile)
  t1 = readImage(ft1file)
  t2 = readImage(ft2file)
  t3 = readImage(ft3file)
  r1 = readImage(fr1file)
  r2 = readImage(fr2file)
  r3 = readImage(fr3file)
  sk = readSkins(fsksubs)
  he = HorizonExtractor()
  rs = array(r1,r2,r3)
  ts = array(t1,t2,t3)
  gx = gain(gx)
  #he.setHorizon(xs,gx)
  #plot3(gx,qg=qg)
  fls = fillfloat(0.0,n1,n2,n3)
  FaultSkin.setValueOnFaults(0.8,sk,fls)
  scs = mul(100,fls)
  gxs = sub(gx,scs)
  slices = [150,194,212]
  qgs = []
  cs = [1.0,0.75,0.5,0.25,0.0]
  for hi,h in enumerate(hs):
    qgi = []
    ws = he.findInUnfaultSpace(h,rs)
    xs = he.findInCurrentSpace(h,ws,ts,rs)
    qg = he.applyForTg(cs[hi],xs,gx,sk)
    '''
    qgi = qgi.append(qg)
    plot3(gxs,g=fls,qg=qgi,cmin=0.01,cmax=1.0,cmap=jetFillExceptMin(1.0),
          slices=slices,cint=0.1,png="horizon"+str(hi))
    '''
    qgs.append(qg)
  plot3(gxs,g=fls,qg=qgs,cmin=0.01,cmax=1.0,cmap=jetFillExceptMin(1.0),
      slices=slices,cint=0.1,png="horizon"+str(hi))

  plot3(gxs,g=fls,cmin=0.01,cmax=1.0,cmap=jetFillExceptMin(1.0),
        slices=slices,cint=0.1,png="seismic")

def goShifts():
  gx = readImage(gxfile)
  gw = readImage(gwfile)
  t1 = readImage(ft1file)
  t2 = readImage(ft2file)
  t3 = readImage(ft3file)
  gs = zerofloat(n1,n2,n3)
  gu = zerofloat(n1,n2,n3)
  rs = array(t1,t2,t3)
  uf = Unfault(6.0,6.0)
  us = uf.convertShifts(20,rs)
  uf.applyShifts(rs,gx,gs)
  uf.applyShifts(us,gx,gu)

  plot3(gw)
  plot3f(gs,clab="use sx")
  plot3f(gu,clab="use ru")
  plot3f(gw,clab="flatten")
  plot3f(gx,a=t1,amin=-10,amax=10,amap=jetFill(0.3),alab="sx1",aint=2.0)
  plot3f(gx,a=us[0],amin=-10,amax=10,amap=jetFill(0.3),alab="ru1",aint=2.0)
  plot3f(gx,a=t2,amin=-2,amax=2,amap=jetFill(0.3),alab="sx2",aint=2.0)
  plot3f(gx,a=us[1],amin=-2,amax=2,amap=jetFill(0.3),alab="ru2",aint=2.0)
  plot3f(gx,a=t3,amin=-2,amax=2,amap=jetFill(0.3),alab="sx3",aint=2.0)
  plot3f(gx,a=us[2],amin=-2,amax=2,amap=jetFill(0.3),alab="ru3",aint=2.0)

def goSlices():
  gx = readImage(gxfile)
  gw = readImage(gwfile)
  fw = readImage(fwsfile)
  t1 = readImage(ft1file)
  t2 = readImage(ft2file)
  t3 = readImage(ft3file)
  r1 = readImage(fr1file)
  r2 = readImage(fr2file)
  r3 = readImage(fr3file)
  fw = smooth(fw)
  gx = gain(gx)
  fw = gain(fw)
  gw = gain(gw)
  sks = readSkins(fsksubs)
  skl = readSkins(fslbase)
  skp = readSkins(fspbase)
  fls = like(gx)
  mark = -100
  fss = fillfloat(mark,n1,n2,n3)
  FaultSkin.getLikelihood(sks,fls)
  FaultSkin.getThrow(mark,skl,fss)
  d1 = 0.002
  mul(t1,d1*1000,t1)
  mul(fss,d1*1000,fss)
  ts = array(t1,t2,t3)
  rs = array(r1,r2,r3)
  he = HorizonExtractor()
  cs = he.computeCompositeShifts(ts,rs)
  c1 = cs[0]
  mul(c1,d1*1000,c1)
  print min(c1)
  print max(c1)
  plot3(gx,skins=skp,smax=20,png="faults")
  plot3f(gx,a=fls,amin=0.01,amax=0.8,
        amap=jetFillExceptMin(1.0),alab="Fault likelihood",aint=0.1,png="fls")
  plot3f(gx,a=fss,amin=-0.01,amax=max(fss)-10,
        amap=jetFillExceptMin(1.0),alab="Fault throw (ms)",
        aint=5.0,png="throws")
  plot3f(fw,a=fss,amin=-0.01,amax=max(fss)-10,
        amap=jetFillExceptMin(1.0),alab="Fault throw (ms)",
        aint=5.0,png="unfault")
  plot3f(gx,a=t1,amin=-24,amax=20,
        amap=jetFill(0.3),alab="Vertical shift (ms)",
        #amap=jetFill(0.3),alab="Vertical component of unfaulting shift (ms)",
        aint=5.0,png="unfaultShifts")
  plot3f(gx,a=c1,amin=-60,amax=90,
        amap=jetFill(0.3),alab="Vertical component of composite shift (ms)",
        aint=20.0,png="unfoldShifts")
  plot3f(fw,clab="Amplitude",aint=0.5,png="unfaultSlices")
  plot3f(gw,clab="Amplitude",aint=0.5,png="unfoldSlices")
  #plot3(gx,c1,cmin=-60,cmax=80,cmap=jetFill(0.3),slices=[93,190,175],
  #      clab="Vertical shift (samples)", cint=20.0,png="compositeShifts3D")


def goSubset():
  gx = readImage(gxfile)
  sks = readSkins(fskgood)
  gx = gain(gx)
  f1,f2,f3=0,  7, 28
  l1,l2,l3=119,206,227
  n1,n2,n3=120,200,200
  gxs = copy(n1,n2,n3,f1,f2,f3,gx)
  writeImage(gxsfile,gxs)
  fss = FaultSkinSub(f1,f2,f3,l1,l2,l3)
  sks = fss.getSubSkin(sks)
  '''
  for iskin,skin in enumerate(sks):
    plot3(gxs,skins=[skin],clab="skin"+str(iskin))
  '''
  plot3s(gxs,skins=sks,png="subSkins")
  sks = [sks[0],sks[1],sks[4],sks[5],sks[6],sks[7],sks[8],sks[9],sks[10],sks[11],
        sks[12],sks[13],sks[14],sks[18],sks[20],sks[22],sks[24]]
  removeAllSkinFiles(fsksubs)
  writeSkins(fsksubs,sks)
  fls = like(gxs)
  FaultSkin.getLikelihood(sks,fls)
  scs = mul(100,fls)
  #gxs = sub(gxs,scs)
  plot3s(gxs,fls,cmin=0.01,cmax=0.8,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",cint=0.2,png="subFl")
  plot3s(gxs,skins=sks,png="subSkins")


def goTeaser():
  gx = readImage(gxfile)
  gx = gain(gx)
  p2 = readImage(p2file)
  p3 = readImage(p3file)
  fl = readImage(flfile)
  fp = readImage(fpfile)
  ft = readImage(ftfile)
  fs = FaultSkinner()
  fs.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
  fs.setMinSkinSize(minSkinSize)
  cells = fs.findCells([fl,fp,ft])
  f1,f2,f3=0,  48,122
  l1,l2,l3=68,147,221
  gxsub = copy(69,100,100,f1,f2,f3,gx)
  fss = FaultSkinSub(f1,f2,f3,l1,l2,l3)
  cellsub = fss.getSubCells(cells)
  sk = readSkins(fskgood)
  skinsub = [sk[0],sk[1],sk[5],sk[6],sk[7],sk[8],sk[10]]
  sk = fss.getSubSkin(skinsub)
  skinsub = [sk[0],sk[1],sk[2],sk[3],sk[4],sk[5],sk[7]]
  cellsub = FaultSkin.getCells(skinsub)
  flsub = like(gxsub)
  FaultSkin.getLikelihood(skinsub,flsub)
  scs = mul(100,flsub)
  gxs = sub(gxsub,scs)
  plot3s(gxs,flsub,cmin=0.01,cmax=0.8,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",cint=0.2,png="subFl")
  plot3s(gxsub,cells=cellsub,png="subCells")
  plot3s(gxsub,skins=skinsub,png="subSkins")

def goReskin():
  useOldCells=True
  gx = readImage(gxfile)
  fl = readImage(flfile)
  sk = readSkins(fskbase)
  if not plotOnly:
    fsx = FaultSkinnerX()
    cells = fsx.resetCells(sk)
    fsx = FaultSkinnerX()
    fsx.setParameters(20.0,10.0,2.0)
    fsx.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsx.setMinSkinSize(minSkinSize)
    fsx.setMaxPlanarDistance(0.5)
    fsx.setSkinning(useOldCells)
    sks = fsx.findSkinsXX(cells,fl)
    removeAllSkinFiles(fskgood)
    writeSkins(fskgood,sks)
  skins = readSkins(fskgood)
  plot3(gx,skins=skins)

def goSmooth():
  print "goSmooth ..."
  gx = readImage(gxfile)
  gx = gain(gx)
  if not plotOnly:
    flstop = 0.01
    fsigma = 8.0
    gx = readImage(gxfile)
    sks = readSkins(fsksubs)
    flt = zerofloat(n1,n2,n3)
    FaultSkin.getLikelihoods(sks,flt)
    sigma1,sigma2,sigma3,pmax = 16.0,2.0,2.0,2.0
    p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,gx)
    writeImage(p2file,p2)
    writeImage(p3file,p3)
    writeImage(epfile,ep)
    gsx = FaultScanner.smooth(flstop,fsigma,p2,p3,flt,gx)
    writeImage(gsxfile,gsx)
  else:
    gsx = readImage(gsxfile)
  plot3(gx,clab="Amplitude",png="gx")
  plot3(gsx,clab="Amplitude",png="gsx")

def goSlip():
  print "goSlip ..."
  gx = readImage(gxfile)
  if not plotOnly:
    skins = readSkins(fsksubs)
    plot3(gx,skins=skins,png="skinsfl")
    gsx = readImage(gsxfile)
    p2 = readImage(p2file)
    p3 = readImage(p3file)
    fsl = FaultSlipper(gsx,p2,p3)
    fsl.setOffset(2.0) # the default is 2.0 samples
    fsl.setZeroSlope(False) # True only to show the error
    fsl.computeDipSlips(skins,minThrow,maxThrow)
    print "  dip slips computed, now reskinning ..."
    print "  number of skins before =",len(skins),
    fsk = FaultSkinner() # as in goSkin
    fsk.setGrowLikelihoods(lowerLikelihood,upperLikelihood)
    fsk.setMinSkinSize(minSkinSize)
    fsk.setMinMaxThrow(-1.0,maxThrow)
    skins = fsk.reskin(skins)
    removeAllSkinFiles(fspbase)
    writeSkins(fspbase,skins)
    plot3(gx,skins=skins,smax=20,png="skinss1")
    fsk.setMinMaxThrow(minThrow,maxThrow)
    skins = fsk.reskin(skins)
    print ", after =",len(skins)
    removeAllSkinFiles(fslbase)
    writeSkins(fslbase,skins)
  else:
    skins = readSkins(fslbase)
  plot3(gx,skins=skins,png="skinsfl")
  plot3(gx,skins=skins,smax=20,png="skinss1")

def unfault():
  smark = -999.999
  gx = readImage(gxfile)
  skins = readSkins(fslbase)
  fsl = FaultSlipper(gx,gx,gx)
  s1,s2,s3=fsl.getDipSlips(skins,smark)
  s1,s2,s3 = fsl.interpolateDipSlips([s1,s2,s3],smark)
  gw = fsl.unfault([s1,s2,s3],gx)
  plot3(gx)
  plot3(gw)

def goUnfaultS():
  if not plotOnly:
    gx = readImage(gxfile)
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(8.0,2.0,2.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    sks = readSkins(fsksubs)
    fls = like(gx)
    FaultSkin.getLikelihood(sks,fls)
    scs = mul(100,fls)
    plot3s(gx,fls,cmin=0.01,cmax=0.8,cmap=jetFillExceptMin(1.0),
        clab="Fault likelihood",cint=0.2,png="subFl")
    plot3s(gx,skins=sks)

    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.screenPoints(wp)

    uf = UnfaultS(4.0,2.0)
    uf.setIters(100)
    uf.setTensors(et)
    mul(sp[3][0],10,sp[3][0])
    [t1,t2,t3] = uf.findShifts(sp,wp)
    plot3(wp)
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwsfile,fw)
    writeImage(ft1file,t1)
    writeImage(ft2file,t2)
    writeImage(ft3file,t3)

  else :
    gx = readImage(gxfile)
    fw = readImage(fwsfile)
    t1 = readImage(ft1file)
    t2 = readImage(ft2file)
    t3 = readImage(ft3file)
  gx = gain(gx)
  fw = gain(fw)
  plot3(gx)
  plot3(fw,clab="unfaulted")
  plot3(gx,t1,cmin=-10,cmax=10,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,t2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,t3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")

def goUnfaultC():
  if not plotOnly:
    gx = readImage(gxfile)
    fw = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(8.0,2.0,2.0)
    et = lof.applyForTensors(gx)
    et.setEigenvalues(0.001,1.0,1.0)

    wp = fillfloat(1.0,n1,n2,n3)
    ws = fillfloat(1.0,n1,n2,n3)
    cp = fillfloat(0.0,n1,n2,n3)

    skins = readSkins(fslbase)
    fsc = FaultSlipConstraints(skins)
    sp = fsc.controlPoints(ws,wp,cp)

    uf = UnfaultC(4.0,4.0)
    uf.setIters(100)
    uf.setTensors(et)
    [t1,t2,t3] = uf.findShifts(sp,wp,ws)
    uf.applyShifts([t1,t2,t3],gx,fw)
    writeImage(fwcfile,fw)
    writeImage(fc1file,t1)
    writeImage(fc2file,t2)
    writeImage(fc3file,t3)

  else :
    gx = readImage(gxfile)
    fw = readImage(fwcfile)
    t1 = readImage(fc1file)
    t2 = readImage(fc2file)
    t3 = readImage(fc3file)
  gx = gain(gx)
  fw = gain(fw)
  plot3(gx)
  plot3(fw,clab="unfaulted")
  plot3(gx,t1,cmin=-10,cmax=10,cmap=jetFill(0.3),
        clab="Vertical shift (samples)",png="gxs1i")
  plot3(gx,t2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="Inline shift (samples)",png="gxs2i")
  plot3(gx,t3,cmin=-1.0,cmax=1.0,cmap=jetFill(0.3),
        clab="Crossline shift (samples)",png="gxs3i")


def goUnfold():
  if not plotOnly:
    fw = readImage(fwsfile)
    gw = zerofloat(n1,n2,n3)
    u1 = zerofloat(n1,n2,n3)
    u2 = zerofloat(n1,n2,n3)
    u3 = zerofloat(n1,n2,n3)
    ep = zerofloat(n1,n2,n3)
    lof = LocalOrientFilter(2.0, 1.0)
    fw = smoothF(fw)
    lof.applyForNormalPlanar(fw,u1,u2,u3,ep)
    pow(ep,8.0,ep)
    p = array(u1, u2, u3, ep)
    flattener = FlattenerRT(4.0,4.0)
    [r1,r2,r3] = flattener.findShifts(p)
    flattener.applyShifts([r1,r2,r3], fw, gw)
    writeImage(gwfile,gw)
    writeImage(fr1file,r1)
    writeImage(fr2file,r2)
    writeImage(fr3file,r3)
  else :
    gw = readImage(gwfile)
    fw = readImage(fwfile)
  gw = gain(gw)
  fw = gain(fw)
  plot3(fw)
  plot3(gw)

def goDisplay():
  gx = readImage(gxfile)
  cp = readImage(cpfile)
  gw = readImage(gwfile)
  ftc = readImage(ftcfile)
  r1 = readImage(r1tfile)
  '''
  r2 = readImage(r2tfile)
  r3 = readImage(r3tfile)
  '''
  ft1 = readImage(ft1file)
  #fs2 = readImage(fs1file)
  #fs3 = readImage(fs1file)
  hmin,hmax,hmap = -1.0,1.0,ColorMap.GRAY
  plot3(ftc,cmin=hmin,cmax=hmax,cmap=hmap,clab="UnfaultC",png="ftc")
  plot3(gw,cmin=hmin,cmax=hmax,cmap=hmap,clab="Unfault",png="gw")
  plot3(gx,cmin=hmin,cmax=hmax,cmap=hmap,clab="Amplitude",png="gx")
  plot3(cp,cmin=hmin,cmax=hmax,cmap=hmap,clab="ControlPoints",png="cp")
  plot3(gx,r1,cmin=-20.0,cmax=10.0,cmap=hueFill(0.3),
        clab="Vertical shift for unfaulting with constraints",png="gxs1i")
  '''
  plot3(gx,r2,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="r2",png="gxs1i")
  plot3(gx,r3,cmin=-2.0,cmax=2.0,cmap=jetFill(0.3),
        clab="r3",png="gxs1i")
  '''
  plot3(gx,ft1,cmin=0.0,cmax=10.0,cmap=jetFill(0.3),
        clab="Vertical shift for unfaulting",png="gxs1i")

def array(x1,x2,x3=None,x4=None):
  if x3 and x4:
    return jarray.array([x1,x2,x3,x4],Class.forName('[[[F'))
  elif x3:
    return jarray.array([x1,x2,x3],Class.forName('[[[F'))
  else:
    return jarray.array([x1,x2],Class.forName('[[[F'))
def like(x):
  n3 = len(x)
  n2 = len(x[0])
  n1 = len(x[0][0])
  return zerofloat(n1,n2,n3)

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def smooth(x):
  fsigma = 2.0
  flstop = 0.9
  sks = readSkins(fsksubs)
  flt = fillfloat(1.0,n1,n2,n3)
  fsc = FaultSlipConstraints(sks)
  fsc.setValuesNearFaults(0.0,flt)
  sigma1,sigma2,sigma3,pmax = 4.0,1.0,1.0,1.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,x)
  return FaultScanner.smooth(flstop,fsigma,p2,p3,flt,x)
def smoothF(x):
  fsigma = 2.0
  flstop = 0.9
  flt = fillfloat(0.0,n1,n2,n3)
  sigma1,sigma2,sigma3,pmax = 4.0,1.0,1.0,1.0
  p2,p3,ep = FaultScanner.slopes(sigma1,sigma2,sigma3,pmax,x)
  return FaultScanner.smooth(flstop,fsigma,p2,p3,flt,x)

#############################################################################
# graphics

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def jetFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.JET,a)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
def hueFill(alpha):
  return ColorMap.getHue(0.0,1.0,alpha)
def hueFillExceptMin(alpha):
  a = fillfloat(alpha,256)
  a[0] = 0.0
  return ColorMap.setAlpha(ColorMap.getHue(0.0,1.0),a)

def addColorBar(frame,clab=None,cint=None):
  cbar = ColorBar(clab)
  if cint:
    cbar.setInterval(cint)
  cbar.setFont(Font("Arial",Font.PLAIN,32)) # size by experimenting
  cbar.setWidthMinimum
  cbar.setBackground(Color.WHITE)
  frame.add(cbar,BorderLayout.EAST)
  return cbar

def convertDips(ft):
  return FaultScanner.convertDips(0.2,ft) # 5:1 vertical exaggeration

def convertStrikes(fp):
  return FaultScanner.convertStrikes(True,-90.0,fp)

def plot3(f,g=None,qg=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
          xyz=None,cells=None,skins=None,smax=0.0,slices=None,
          links=False,curve=False,trace=False,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  #sf = SimpleFrame(AxesOrientation.XRIGHT_YOUT_ZDOWN)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-2.0,2.0)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-2.0,2.0)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if xyz:
    pg = PointGroup(0.2,xyz)
    ss = StateSet()
    cs = ColorState()
    cs.setColor(Color.YELLOW)
    ss.add(cs)
    pg.setStates(ss)
    #ss = StateSet()
    #ps = PointState()
    #ps.setSize(5.0)
    #ss.add(ps)
    #pg.setStates(ss)
    sf.world.addChild(pg)
  if qg:
    ss = StateSet()
    lms = LightModelState()
    #lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    for qgi in qg:
      qgi.setStates(ss)
      sf.world.addChild(qgi)
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    cmap = ColorMap(0.0,0.8,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,True)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    if not smax:
      ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.5 
    for skin in skins:
      if smax>0.0: # show fault throws
        cmap = ColorMap(-1.0,smax,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForThrow(size,cmap,True)
      else: # show fault likelihood
        cmap = ColorMap(0.0,0.8,ColorMap.JET)
        xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,True)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if curve or trace:
        cell = skin.getCellNearestCentroid()
        if curve:
          xyz = cell.getFaultCurveXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
        if trace:
          xyz = cell.getFaultTraceXyz()
          pg = PointGroup(0.5,xyz)
          sg.addChild(pg)
      if links:
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  if slices:
    k1,k2,k3 = slices
  else:
    #k1,k2,k3 = (370,105,34) # most plots use these
    #k1,k2,k3 = (370,0,0) # most plots use these
    #k1,k2,k3 = (200,208,228) # most plots use these
    k1,k2,k3 = (150,194,212) # most plots use these
    #k1,k2,k3 = (370,150,0) # most plots use these
    k1,k2,k3 = (120,100,100) # most plots use these
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(985,700) # for sch data
    #sf.setSize(837,700) # for fake data
  else:
    sf.setSize(848,700) # for sch data
    #sf.setSize(1048,900) # for sch data
    #sf.setSize(700,700) # for fake data
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setEyeToScreenDistance(3018.87) # for consistency with brooks
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(25.0,20.0)
  ov.setAzimuthAndElevation(-35.0,50.0)
  #ov.setAzimuthAndElevation(150.0,15.0)
  #ov.setAzimuthAndElevation(160.0,65.0)
  ov.setScale(1.3)
  #ov.setTranslate(Vector3(-0.182,-0.238,-0.012))
  #ov.setTranslate(Vector3(-0.190,-0.168,-0.006))
  ov.setTranslate(Vector3(-0.290,-0.308,-0.00))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3s(f,g=None,cmin=None,cmax=None,cmap=None,clab=None,cint=None,
           cells=None,skins=None,links=False,png=None):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  s1,s2,s3=Sampling(n1),Sampling(n2),Sampling(n3)
  sf = SimpleFrame(AxesOrientation.XRIGHT_YIN_ZDOWN)
  cbar = None
  if g==None:
    ipg = sf.addImagePanels(s1,s2,s3,f)
    if cmap!=None:
      ipg.setColorModel(cmap)
    if cmin!=None and cmax!=None:
      ipg.setClips(cmin,cmax)
    else:
      ipg.setClips(-1.5,1.5)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMapListener(cbar)
  else:
    ipg = ImagePanelGroup2(s1,s2,s3,f,g)
    ipg.setClips1(-1.5,1.5)
    if cmin!=None and cmax!=None:
      ipg.setClips2(cmin,cmax)
    if cmap==None:
      cmap = jetFill(0.8)
    ipg.setColorModel2(cmap)
    if clab:
      cbar = addColorBar(sf,clab,cint)
      ipg.addColorMap2Listener(cbar)
    sf.world.addChild(ipg)
  if cbar:
    cbar.setWidthMinimum(120)
  if cells:
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    cmap = ColorMap(0.0,0.8,ColorMap.JET)
    xyz,uvw,rgb = FaultCell.getXyzUvwRgbForLikelihood(0.5,cmap,cells,True)
    qg = QuadGroup(xyz,uvw,rgb)
    qg.setStates(ss)
    sf.world.addChild(qg)
  if skins:
    sg = Group()
    ss = StateSet()
    lms = LightModelState()
    lms.setTwoSide(True)
    ss.add(lms)
    ms = MaterialState()
    ms.setSpecular(Color.GRAY)
    ms.setShininess(100.0)
    ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE)
    ms.setEmissiveBack(Color(0.0,0.0,0.5))
    ss.add(ms)
    sg.setStates(ss)
    size = 2.0
    if links:
      size = 0.5 
    ik = 0
    for skin in skins:
      cmap = ColorMap(0.0,0.8,ColorMap.JET)
      xyz,uvw,rgb = skin.getCellXyzUvwRgbForLikelihood(size,cmap,True)
      qg = QuadGroup(xyz,uvw,rgb)
      qg.setStates(None)
      sg.addChild(qg)
      if links:
        xyz = skin.getCellLinksXyz()
        lg = LineGroup(xyz)
        sg.addChild(lg)
    sf.world.addChild(sg)
  #k1,k2,k3 = (n1,n2,n3) # most plots use these
  k1,k2,k3 = (n1,87,90) # most plots use these
  ipg.setSlices(k1,k2,k3)
  if cbar:
    sf.setSize(985,700) # for sch data
    #sf.setSize(837,700) # for fake data
  else:
    sf.setSize(848,700) # for sch data
  vc = sf.getViewCanvas()
  vc.setBackground(Color.WHITE)
  radius = 0.5*sqrt(n1*n1+n2*n2+n3*n3)
  ov = sf.getOrbitView()
  ov.setEyeToScreenDistance(3018.87) # for consistency with brooks
  ov.setWorldSphere(BoundingSphere(0.5*n1,0.5*n2,0.5*n3,radius))
  #ov.setAzimuthAndElevation(25.0,20.0)
  ov.setAzimuthAndElevation(-45.0,40.0)
  #ov.setAzimuthAndElevation(150.0,15.0)
  #ov.setAzimuthAndElevation(160.0,65.0)
  ov.setScale(1.2)
  #ov.setTranslate(Vector3(-0.182,-0.238,-0.012))
  ov.setTranslate(Vector3(-0.190,-0.168,-0.006))
  sf.setVisible(True)
  if png and pngDir:
    sf.paintToFile(pngDir+png+".png")
    if cbar:
      cbar.paintToPng(137,1,pngDir+png+"cbar.png")

def plot3f(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           clab=None,png=None):
  background = Color.WHITE
  n3 = len(g)
  n2 = len(g[0])
  n1 = len(g[0][0])
  d1,d2,d3 = 0.002,0.008,0.008 # (s,km,km)
  f1,f2,f3 = 0.300,0.0,0.0 # (s,km,km)
  s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  #k1,k2,k3=78,193,175
  #k1,k2,k3=78,193,168
  #k1,k2,k3=63,189,174
  #k1,k2,k3=63,180,166
  k1,k2,k3=63,182,166
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  #pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1, 85)
  pp.setClips(-1.5,1.5)
  if a:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar(clab)
    cb.setInterval(0.5)
  pp.setInterval1(0.1)
  pp.setInterval2(0.3)
  pp.setInterval3(0.3)
  if a:
    pv12 = PixelsView(s1,s2,slice12(k3,a))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,a))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1,a))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(amap)
      if amin!=amax:
        pv.setClips(amin,amax)
    pp.getPixelsView12.Tile.addTiledView(pv12)
    pp.getPixelsView13.Tile.addTiledView(pv13)
    pp.getPixelsView23.Tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(100)
  #pf.setFontSize(18) #for print
  pf.setFontSize(30) #for slices
  #pf.setFontSizeForPrint(1.0,0.8)
  pf.setSize(1150,800)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+png
    pf.paintToPng(360,7.0,png+".png")
  
def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s
 


#############################################################################
run(main)
