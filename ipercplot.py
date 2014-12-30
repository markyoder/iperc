import math
import random
import operator
import matplotlib.pyplot as plt
import ipercPy as ipc
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import scipy
import scipy.stats
import scipy.optimize as spo
import datetime as dtm
import matplotlib.dates as mpd
import pytz
import os

import ANSStools as atp
import eqcatalog as eqp
import linefit as lfp

plt.ion()
#
# some work...
def iperkSizeSpect(nits=1000, width=1024, geom=0, dla=1, fnum=0):
	# spectrum of cluster sizes for a range of lattice sizes.
	ks=[]
	csize=ipc.ipercObj(width, geom, dla)-1
	for i in xrange(nits):
		#csize=ipc.getClustsize();
		ipc.initializeLattice()
		ipc.initializeCluster()
		ipc.makePercolatingCluster()
		csize=ipc.getClustsize()
		ks+=[csize]
		ipc.randSeedPlus()
		#print "size: %d" % csize
		#reload(ipc)
		#
		'''
		plt.figure(fnum)
		plt.clf()
		XY = getXYfromClust(ipc.getClusterSites(), width	)
		plt.plot(XY[0], XY[1], '.')
		fnum+=1
		'''
	#
	ks.sort()
	ns=range(len(ks))
	ns.reverse()
	plt.figure(fnum)
	plt.clf()
	plt.ion()
	plt.plot(ks, ns, '-')
	#
	return ks

def getNNpairs(zs, zscale=None, L=1024, geom=0, dla=1, fnum=0):
	if zs==None:
		zs=getNormalZs(L=L, geom=geom, dla=dla)
	if zscale==None: zscale=float(sum(zs))
	#
	zs=scipy.array(zs)
	dzs=zs[1:]-zs[0:-1]
	#rsqr = 
	#
	
def RBclusters(zs=None, L=1024, geom=0, dla=1, fnum=0, xscale='linear', yscale='linear'):
	if zs==None:
		zs=getNormalZs(L=L, geom=geom, dla=dla)
	if xscale=='lin': xscale='linear'
	if yscale=='lin': yscale='linear'
	zs=scipy.array(zs)
	dzs=zs[1:]-zs[0:-1]
	csize=float(len(zs))
	print "cluster size: %d" % csize	
	#
	# first, what is our time-series?
	f=plt.figure(fnum)
	plt.clf()
	#sumr=[0]
	#for z in zs:
	#	sumr+=[sumr[-1]+z]
	delta_ts=getIntsfromZs(zs)
	sumts=[delta_ts[0]]
	sumzs=[zs[0]]
	for i in xrange(1, len(zs)):
		sumts+=[sumts[-1]+delta_ts[i]]
		sumzs+=[sumzs[-1]+zs[i]]
	#
	#print "sumrLens: %f, %f" % (sumts[-1], len(zs))
	#plt.plot(sumr[1:], zs, ',')
	axdn=f.add_axes([.1, .05, .85, .4])
	axup=f.add_axes([.1, .5, .85, .45], sharex=axdn)
	axup.set_yscale('log')
	axdn.plot(sumts, zs, ',')
	axdn.set_ylabel('$z$')
	axdn.set_xlabel('Weibull time $t_w$')
	axup.plot(sumts, 1.0/scipy.array(delta_ts), ',')
	#ax.plot(sumzs)
	plt.title('Weibull $z$ and interval $\\Delta t_w$ plots')
	#plt.xlabel('Weibull time $t_w$')
	axup.set_ylabel('Weibull interval $\\Delta t_w$')
	#
	# a few definitions.
	# 1) a record-breaking cluster: all record-breaking largest events between recordbreaking smallest events
	# (or vice versa)
	#
	rbclusters_big=[[zs[0]]]
	X=[sumts[0]]
	zbig=zs[0]
	zsmall=zs[0]
	nrb_big=1
	nrb_small=1
	for i in xrange(1,int(csize)):
		z=zs[i]
		# do we start a new cluster?
		if z<zsmall:
			# new cluster.
			zbig=z
			zsmall=z
			rbclusters_big+=[[]]
		# do we include?
		if z>=zbig:
			# new rb biggest.
			rbclusters_big[-1]+=[z]
			zbig=z
			X+=[sumts[i]]
			#nrb_big+=1
		#
	#
	# plot the full sequence.
	allzs=[]
	# there's a better way, but i don't know it off the top of my head...
	for rw in rbclusters_big:
		for elem in rw:
			allzs+=[elem]
		#
	#
	rbclusters_small=[[zs[0]]]
	X2=[sumts[0]]
	zbig=zs[0]
	zsmall=zs[0]
	nrb_big=1
	nrb_small=1
	for i in xrange(1,int(csize)):
		z=zs[i]
		# do we start a new cluster?
		if z>zbig:
			# new cluster.
			zbig=z
			zsmall=z
			rbclusters_small+=[[]]
		# do we include?
		if z<=zsmall:
			# new rb biggest.
			rbclusters_small[-1]+=[z]
			zsmall=z
			X2+=[sumts[i]]
			#nrb_big+=1
		#
	#
	# plot the full sequence.
	allzs=[]
	allzs2=[]
	# there's a better way, but i don't know it off the top of my head...
	for rw in rbclusters_big:
		for elem in rw:
			allzs+=[elem]
	for rw in rbclusters_small:
		for elem in rw:
			allzs2+=[elem]	
	axdn.plot(X, allzs, '.-')
	axdn.plot(X2, allzs2, '.-')
	
def getNTWzs(fname='synthcats/rundleNTW.dat'):
	# get NTW data (Rundle et al. 2011,2012). convert these into an i-perc catalog of z-values (strengths) for comparison to
	# synthetic data.
	f=open(fname)
	t=[]
	zs=[]
	#
	# this file is a mess, but parsable...
	#i=0
	for rw in f:
		#rw_len=len(rw)
		#if i>20: break
		#i+=1
		if rw[0] in ('#', '\n', '\t'): continue
		rws=map(float, rw.split())
		#print len(rws), " ** ", rws	
		#if len(rws)!=6: continue
		#print rws
		
		#
		if rws[0]<=100:
			# it's a prob.
			zs+=rws
		if rws[0]>100:
			t+=rws
	#
	f.close()
	zs=scipy.array(zs)/100.
	#
	return [t, zs]
	
def NTWset(fname='synthcats/rundleNTW.dat'):
	#
	tz=getNTWzs(fname='synthcats/rundleNTW.dat')
	t=tz[0]
	zs=tz[1]
	#
	plt.figure(0)
	plt.clf()
	plt.plot(t, zs, '-')
	#
	dzs=zs[1:]-zs[0:-1]
	#
	plt.figure(1)
	plt.clf()
	plt.hist(dzs, bins=200, log=False, normed=True)
	
	
	#
	return [t, zs]

def plotDeltaRset(zs=None, L=1024, geom=0, dla=1, fnum=0, xscale='linear', yscale='linear'):
	if zs==None:
		zs=getNormalZs(L=L, geom=geom, dla=dla)
	if xscale=='lin': xscale='linear'
	if yscale=='lin': yscale='linear'
	zs=scipy.array(zs)
	dzs=zs[1:]-zs[0:-1]
	csize=float(len(zs))
	print "cluster size: %d" % csize
	# first, just plot R(t).
	'''
	plt.figure(fnum)
	plt.clf()
	plt.plot(zs, ',')
	plt.xlabel('$t_{natural}$')
	plt.ylabel('$z$')
	#
	# now, a Poisson plot:
	plt.figure(fnum+1)
	plt.clf()
	plt.plot(zs[0:-1], zs[1:], ',')
	plt.xlabel('$r_{i-1}$')
	plt.ylabel('$z_i$')
	#
	# now, plot r_i-r_{i-1}:
	plt.figure(fnum+2)
	plt.clf()
	plt.plot(dzs, ',')
	plt.ylabel('$\\Delta z_i$ $(r_i - z_{i-1} )$')
	plt.xlabel('$t_{natural}$')
	#
	# now, plot [z_i-r_{i-1}][z_i]:
	plt.figure(fnum+3)
	plt.clf()
	plt.plot(zs[1:], dzs, ',')
	plt.ylabel('$\\Delta z_i$ $(r_i - z_{i-1} )$')
	plt.xlabel('$z_i$')
	#
	'''
	#
	# and now, the histogram:
	plt.figure(11)
	plt.clf()
	plt.hist(dzs, bins=400, normed=True, log=False)
	plt.figure(fnum+4)
	plt.clf()
	#plt.hist(dzs, bins=10, normed=True)
	hlogged=False
	plt.hist(dzs, bins=25, normed=True, log=hlogged)
	plt.hist(dzs, bins=50, normed=True, log=hlogged)
	plt.hist(dzs, bins=100, normed=True, log=hlogged)
	plt.hist(dzs, bins=250, normed=True, log=hlogged)
	nhrh=200
	a=plt.hist(dzs, bins=nhrh, normed=True, log=hlogged, histtype='step')
	hireshist=a
	#
	#random:
	xr1=((scipy.arange(100)-100)/100.)
	#xr1.reverse()
	xr2=((scipy.arange(101))/100.)
	yr1=xr1+1.0
	yr2=xr2.tolist()
	yr2.reverse()
	#plt.plot(xr1 + xr2, yr1+yr2, '-')
	plt.plot(xr1, yr1, 'r', alpha=11)
	plt.plot(xr2, yr2, 'r', alpha=11)
	#
	a=scipy.array(a)
	hwidth1=float(a[1][1]-a[1][0])
	#thisPs=a[0]*hwidth1
	#thisDrs=a[1][1:]
	print "sum Probs[all, left, right] (500): %f, %f, %f" % (sum(a[0])*hwidth1, sum(a[0][0:(len(a[0])/2)])*hwidth1, sum(a[0][-(len(a[0])/2):])*hwidth1)
	#
	f=plt.figure(fnum+5)
	plt.clf()
	#sumr=[0]
	#for z in zs:
	#	sumr+=[sumr[-1]+z]
	delta_ts=getIntsfromZs(zs)
	sumts=[delta_ts[0]]
	sumzs=[zs[0]]
	for i in xrange(1, len(zs)):
		sumts+=[sumts[-1]+delta_ts[i]]
		sumzs+=[sumzs[-1]+zs[i]]
	#
	#print "sumrLens: %f, %f" % (sumts[-1], len(zs))
	#plt.plot(sumr[1:], zs, ',')
	axdn=f.add_axes([.1, .05, .85, .4])
	axup=f.add_axes([.1, .5, .85, .45], sharex=axdn)
	axup.set_yscale('log')
	axdn.plot(sumts, zs, ',')
	axdn.set_ylabel('$z$')
	axdn.set_xlabel('Weibull time $t_w$')
	axup.plot(sumts, 1.0/scipy.array(delta_ts), ',')
	#ax.plot(sumzs)
	plt.title('Weibull $z$ and interval $\\Delta t_w$ plots')
	#plt.xlabel('Weibull time $t_w$')
	axup.set_ylabel('Weibull interval $\\Delta t_w$')
	
	# let's look at delta-r >0 and <0 separately:
	'''
	dzgt=[]
	dzlt=[]
	#
	for dz in dzs:
		if dz>=0: dzgt+=[dz]
		if dz<0: dzlt+=[dz]
	print "num dz<0: %f, num dz>=0: %f, (%f)" % (len(dzlt), len(dzgt), float(len(dzgt))/float(len(dzlt)))
	print "sum dz<0: %f, sum dz>=0: %f, (%f)" % (sum(dzlt), sum(dzgt), float(sum(dzgt))/float(sum(dzlt)))
	dzgt=scipy.array(dzgt)
	dzlt=-scipy.array(dzlt)
	'''
	wdth=hireshist[1][1]-float(hireshist[1][0])
	hgt,hlt,xgt,xlt = [], [], [], []
	for i in xrange(len(hireshist[0])):
		x=hireshist[1][i]
		#
		if x>(wdth):
			hgt+=[hireshist[0][i]]
			xgt+=[x]
		if x<(-wdth):
			hlt+=[hireshist[0][i]]
			xlt+=[-x]
			#
		#
	#hgt = hireshist[0][(nhrh/2):]
	#hlt = hireshist[0][0:(nhrh/2)]
	#xgt = hireshist[1][(nhrh/2):-1] + wdth/2.
	#xlt = -(hireshist[1][0:(nhrh/2)] + wdth/2.)
	#xlt.reverse()
	#
	plt.figure(fnum+6)
	plt.clf()
	#h1=plt.hist(dzgt, bins=100, histtype='step', normed=True)
	#h2=plt.hist(dzlt, bins=100, histtype='step', normed=True)
	plt.plot(xgt, hgt, '-')
	plt.plot(xlt, hlt, '-')
	print "+/- balance (approximate): %f/%f" % (sum(hgt), sum(hlt))
	#
	#print "sumh1: %f, sum2: %f" % (sum(h1[0]*float(h1[1][1]-h1[1][0])), sum(h2[0])*float(h2[1][1]-h2[1][0]))
	hwidth6=wdth
	#
	########################
	plt.figure(fnum+7)
	plt.clf()
	ax=plt.gca()
	ax.set_yscale(yscale)
	ax.set_xscale(xscale)
	ax2=ax.twinx()
	ax2.set_yscale(yscale)
	ax2.set_xscale(xscale)
	#plt.plot(h1[1][1:], h1[0]*hwidth6, 'b.--', label='$P>0$')
	#plt.plot(h2[1][1:], h2[0]*hwidth6, 'g.--', label='$P<0$')
	ax.plot(xgt, hgt, 'bo--', label='$P>0$')
	ax.plot(xlt, hlt, 'go--', label='$P<0$')
	#
	# Linear Fit:
	#lf1=lfp.linefit([h2[1][1:], h2[0]])
	lf1 = lfp.linefit([xlt,hlt])
	plsq=lf1.doFit()[0]
	print "lin fit to left: ", plsq, lf1.meanVar()
	ax.plot(xlt, (plsq[0] + plsq[1]*scipy.array(xlt)), 'g-', label='$P<$ linfit')
	#
	# hard fits:
	#ly1=map(math.log, h1[0]*hwidth6)
	#lx1=map(math.log, h1[1][1:])
	lx1, ly1, X, Y = [],[], [],[]
	# we need some manacured data (basically remove zeros).
	#for i in xrange(len(h1[0])):
	for i in xrange(len(hgt)):
		if hgt[i]<=0.0 or xgt[i]<=0.0: continue
		Y+=[hgt[i]]
		X+=[xgt[i]]
	ly1 = scipy.array(map(math.log10, Y))
	lx1 = scipy.array(map(math.log10, X))
	#
	Y=scipy.array(Y)
	X=scipy.array(X)
	#	
	lf2=lfp.linefit([X, ly1])
	lf3=lfp.linefit([lx1,ly1])
	#
	plsq2=lf2.doFit()[0]	# exponential
	plsq3=lf3.doFit()[0]	# PL
	plsqPLexp=spo.leastsq(plexpres, scipy.array([-2., -1., -1.]), args=(Y, X, None) )[0] 
	plsqPLtru=spo.leastsq(pltruncres, scipy.array([-2., .25, 1.]), args=(Y, X, .5) )[0]			# with x0=.5 is fixed
	plsqPLtru2=spo.leastsq(pltruncres, scipy.array([-2., .25, 1.]), args=(Y, X, None) )[0]
	#
	# remember, pltruncres, etc. return log values, particularly log(x0), so log(.5)~-.3
	#
	print "plsqExp:", plsq2
	print "plsqPLa: ", plsq3
	print "plsqPLexp: ", plsqPLexp
	print "plsqPLtru: ", plsqPLtru
	print "plsqPLtru2: ", plsqPLtru2
	#
	#print "plsqPlexp: ", plsqPLexp
	#
	#print "log-lin fit to right: ", plsq2[0], lf2.meanVar()
	#print "log-log fit to right: ", plsq3, len(plsq3)
	#
	beta=scipy.ones(len(X))*plsq3[1]
	ypl = (10**plsq3[0])*(X**plsq3[1])
	yexp = 10**(plsq2[0] + plsq[1]*X)
	#
	yplexp = plexp(plsqPLexp, X)
	ypltru = pltrunc(plsqPLtru, X)
	ypltru2 = pltrunc(plsqPLtru2, X)
	deriv = (10**plsqPLtru[0])*(plsqPLtru[1])*X**(-(1.0 + plsqPLtru[1]))
	
	#plt.plot(X, map(math.exp, (plsq2[0] + plsq2[1]*X)), 'b-')
	ax.plot(X, ypl, 'r.-.', label='PLfit')
	ax.plot(X, yexp, 'c.-.', label='exp fit')
	ax.plot(X, yplexp, 'y.-.', label='plexp fit')
	ax.plot(X, ypltru, 'k.-.', alpha=.6, label='pltrunc fit')
	ax.plot(X, ypltru2, 'b.-.', alpha=.6, label='pltrunc fit2')
	ax2.plot(X, deriv, 'm--', alpha=.7, label='deriv(PLTrunc)')
	ax.legend(loc='best')
	ax2.legend(loc='best')
	#
	return zs

def plres(prams, y,x):
	A=10.0**prams[0]
	beta = scipy.ones(len(x))*prams[1]
	#
	ys=A*scipy.array(map(pow, x, beta))
	#	
	lyt = map(math.log10, ys)
	return (scipy.array(lyt)-scipy.array(map(math.log10, y)))

def pltrunc(prams, x, x0=None):
	A=10**prams[0]
	if x0!=None:
		#prams[2]=math.log10(.5)
		prams[2]=math.log10(x0)
	beta = float(prams[1])
	lx0 = float(prams[2])
	#beta = scipy.ones(len(x))*prams[1]
	#x0=scipy.ones(len(x)) * prams[2]
	#x=scipy.array(x)
	#
	#return scipy.array(A*(x**(-beta) - x0**(-beta)))
	#return scipy.array(A*(x**(-prams[1]) - 10**(-prams[2]*prams[1])))
	return scipy.array(A*(x**(-beta) - 10**(-beta*lx0)))
	
def pltruncres(prams, y,x, x0=None):
	yth=pltrunc(prams, x, x0)
	ly = []
	lyt=[]
	#print y, yth
	for i in xrange(len(y)):
		if y[i]<=0. or yth[i]<=0.: continue
		ly+=[math.log10(y[i])]
		lyt+=[math.log10(yth[i])]
	ly=scipy.array(ly)
	lyt=scipy.array(lyt)
	#
	#return (map(math.log10, yth) - scipy.array(map(math.log10, y)))
	#print "trlens: ", len(y), len(ly), len(lyt)
	return ly-lyt

def plexp(prams, x, x00=None):
	A=10.0**prams[0]
	if x00!=None:
		#prams[2]=math.log10(.5)
		prams[2]=math.log10(x00)
	#prams[2]=math.log10(.5)
	beta = prams[1]
	x0 = 10.**float(prams[2])
	x=scipy.array(x)
	#
	ypl=x**(-beta)
	yexp=scipy.array(map(math.exp, -x/x0))
	
	return A*ypl*yexp
	
def plexpres(prams, y,x, x0=None):
	#
	ymodel=plexp(prams, x)
	#
	return scipy.array(map(math.log10, ymodel)) - scipy.array(map(math.log10,y))

def pow10(x):
	return 10**x
def powe(x):
	return (math.e)**x	
def powalpha(X, alpha):
	Y=[]
	for x in X:
		Y+=[x**alpha]
	return scipy.array(Y)
def plotZbursts(L=1024, geom=0, dla=1, fnum=0, thresh=.35, nits=10):
	plt.figure(fnum)
	plt.clf()
	
#def getNormalZs(zs=None, L=1024, geom=0, dla=1, thresh=.35):
def getNormalZs(zs=None, L=1024, geom=0, dla=1):
	# this assumes zs are raw from the lattice. nominally the proper normalization
	# would be to divide by the max value or some known factor thereof.
	if zs==None:
		a=ipc.init(L, geom)
		ipc.setdlamode(dla)
		#nnn=ipc.getnNN()
		a = ipc.initializeLattice()
		b = ipc.initializeCluster()
		sz=ipc.makePercolatingCluster(dla)
		cbonds=ipc.getClusterBonds()
		zs=map(operator.itemgetter(2), cbonds)
	nnn=ipc.getnNN()
	zs = scipy.array(zs)
	zs = 1.0 - (zs-10.0)/(float(L)*float(L)*float(nnn))
	return zs
	
def getZbursts(zs=None, L=1024, geom=0, dla=1, fnum=0, thresh=.35):
	if zs==None:
		zs=getNormalZs(L=L, geom=geom, dla=dla)
		'''
		a=ipc.init(L, geom)
		ipc.setdlamode(dla)
		nnn=ipc.getnNN()
		a = ipc.initializeLattice()
		b = ipc.initializeCluster()
		sz=ipc.makePercolatingCluster(dla)
		cbonds=ipc.getClusterBonds()
		zs=map(operator.itemgetter(2), cbonds)
		zs = scipy.array(zs)
		zs = 1.0 - (zs-10.0)/(float(L)*float(L)*float(nnn))
		'''
	#
	bursts=[0]
	
	#
	for i in xrange(len(zs)):
		thisz=zs[i]
		if thisz>=thresh: bursts[-1]+=1
		#
		if thisz<thresh and bursts[-1]>0: bursts+=[0]
	#
	bursts.sort()
	Ns=range(1, len(bursts)+1)
	#
	plt.figure(0)
	plt.clf()
	plt.plot(zs, '.-.', alpha=.5)
	#
	plt.figure(1)
	plt.clf()
	#plt.plot(bursts, Ns, '.-')
	plt.hist(bursts, bins=25, histtype='bar', log=True)
	#
	return bursts
	

def displZs(L=1024, geom=0, dla=1, fnum=0, nits=1, nskips=0, nbins=100):
	a = ipc.init(L,geom)
	ipc.setdlamode(dla)
	#a = ipc.initializeLattice()
	nnn=ipc.getnNN()
	if nskips==None: nskips=L*5
	plt.figure(fnum)
	plt.clf()
	#plt.figure(fnum+1)
	#plt.clf()
	for boogy in xrange(nits):
		a = ipc.initializeLattice()
		print "lattice initialized..."
		b = ipc.initializeCluster()
		print "lat-len: %d, clust-len: %d" % (ipc.getLatSize(), len(ipc.getClusterSites()))
		#ipc.setdlamode(dla)
		#
		sz=ipc.makePercolatingCluster(dla)
		#zs = ipc.getDisplacedZs()		# these will just be random, since we're doing a bond-method at this point.
		cbonds=ipc.getClusterBonds()
		zs=map(operator.itemgetter(2), cbonds)
		zs = scipy.array(zs)
		zs = 1.0 - (zs-10.0)/(float(L)*float(L)*float(nnn))
		#zs = (zs-10.0)/(float(L)*float(L)*float(nnn))
		#
		plt.figure(fnum)
		plt.hist(zs[nskips:], nbins, histtype='step', normed=True)
		#
		# see justPlotDispZs() for the following plot:
		#plt.figure(fnum+1)
		#zs2=scipy.array(zs[nskips:]).copy()
		#zs2.sort()
		#plt.plot(zs2[1:], scipy.array(zs2[1:])-scipy.array(zs2[:-1]), '.-.')
		#
		print "size: %d" % sz
	#
	return zs

def justPlotDispZs(zs, fnum=0, nskips=0):
	#if nskips==None: nskips=L*5
	# this can be adapted to produce histograms (trivially). right now, it's
	# experimenting with an approach analogous to intervals-no-binning for omori analysis.
	# basically, we're looking at \Delta r (r) or (1/Delta r)(r)
	plt.figure(fnum)
	plt.clf()
	zs2=scipy.array(zs[nskips:]).copy()
	zs2.sort()
	plt.plot(zs2[1:], (scipy.array(zs2[1:])-scipy.array(zs2[:-1])), '.')

def ipercCorrLen(nits=5, width=1024, geom=0, dla=1, fnum=0, nbins=50, donorm=False, i0=0, dt0=10, D=1.82):
	ipc.init(width, geom)
	ipc.setdlamode(dla)
	Cs=[]
	#plt.figure(fnum)
	#plt.clf()
	#plt.figure(fnum+1)
	#plt.clf()
	rseed = ipc.getrandSeed()
	
	for i in xrange(fnum, fnum+5):
		plt.figure(i)
		plt.clf()
	
	drs1Total=[]
	drs2Total=[]
	
	for i in xrange(nits):
		#print "i=", i
		# &thislatlen, &latgeom, &dlamode, &randseed
		#ipc.globLat = ipc.ipercObj(width, geom, dla, rseed)
		rseed+=1
		#
		#
		ipc.init(width, geom, rseed)
		ipc.initializeLattice()
		ipc.initializeCluster()
		ipc.setdlamode(dla)
		#
		print "Initial cluster: ", ipc.getClusterSites(), "[%d]" %i, " [%d]" % ipc.getrandSeed()
		#print "now, percolate..."
		sz=ipc.makePercolatingCluster(dla)
		#
		Cmap=ipc.getClustermap()
		drs1 = plotDeltaRs(Cmap, i0=i0, dt=1, fnum=fnum+1, doclf=False, symb='-', label='dt=1, %d' % i)
		drs2 = plotDeltaRs(Cmap, i0=i0, dt=dt0, fnum=fnum+1, doclf=False, symb='--', label='dt=%d, %d' % (dt0, i))
		#
		drs1Total+=drs1
		drs2Total+=drs2
		#
		plt.figure(fnum+2)
		plt.title('$\\Delta r(i+1, i)$ distributions$')
		#plt.clf()
		ax=plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		h11=plt.hist(drs1, bins=nbins, histtype='step', cumulative=False, color='b', normed=donorm)
		h12=plt.hist(drs2, bins=nbins, histtype='step', cumulative=False, color='g', normed=donorm)
		#
		# normalized scaling function(s):
		X11, Y11, X12, Y12 = [], [],[],[]
		#
		for i in xrange(len(h11[0])):
			Y11+=[h11[0][i]*h11[1][i]]	# r*N
			X11+=[(h11[1][i]**D)/1.0]
		for i in xrange(len(h12[0])):
			Y12+=[h12[0][i]*h12[1][i]]	# r*N
			X12+=[(float(h12[1][i])**D)/float(dt0)]
		plt.figure(fnum+4)
		plt.title("scaling functions")
		#plt.clf()
		plt.loglog(X11,Y11, '-')
		plt.loglog(X12, Y12, '-')
		
		#plt.figure(fnum+3)
		#ax=plt.gca()
		#ax.set_xscale('log')
		#ax.set_yscale('log')
		#plt.title('Histograms of drs')
		#h21=plt.hist(drs1, bins=nbins, histtype='step', cumulative=True, color='b', normed=donorm)
		#h22=plt.hist(drs2, bins=nbins, histtype='step', cumulative=True, color='g', normed=donorm)
		#
		#fnum+=1
		ipc.randSeedPlus()
	#
	plt.figure(fnum)
	plt.clf()
	Z = ipc.getClusterBonds()
	z = plotClustBonds(Z, fnum=(fnum+0), w=width, geom=geom, dla=dla)
	#xyzv = ipc.getClustermap()
	#x=map(operator.itemgetter(1), xyzv)
	#y=map(operator.itemgetter(2), xyzv)
	#plt.plot(x,y,'.')
	#
	plt.figure(5)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	h1t=plt.hist(drs1Total, bins=nbins, histtype='step', cumulative=False)
	h2t=plt.hist(drs2Total, bins=nbins, histtype='step', cumulative=False)
	#
	X11, Y11, X12, Y12 = [], [],[],[]
	#D=1.82
	for i in xrange(len(h11[0])):
		Y11+=[h1t[0][i]*h1t[1][i]]	# r*N
		X11+=[(h1t[1][i]**D)/1.0]
	for i in xrange(len(h12[0])):
		Y12+=[h2t[0][i]*h2t[1][i]]	# r*N
		X12+=[(float(h2t[1][i])**D)/float(dt0)]
	plt.figure(6)
	plt.clf()
	plt.title("Collected scaling functions")
	plt.loglog(X11,Y11, '-')
	plt.loglog(X12, Y12, '-')
		
	plt.figure(fnum+1)
	plt.title('plotDeltaRs()')
	plt.legend(loc='best')
	#
	#plt.figure(fnum+2)
	#h1=plt.hist(drs, bins=10)
	#	
	return drs1
	
def plotDeltaRs(cmap=None, i0=0, dt=1, fnum=0, doclf=True, symb='.-', label=''):
	# CDF of corr-len
	'''
	if zs==None: zs=getNormalZs(L=L, geom=geom, dla=dla)
	zs=zs[i0:]
	cmap=ipc.getClustermap()
	'''
	if cmap==None:
		zs=getNormalZs()
		cmap=ipc.getClustermap()
	
	drs = getDeltaRs(cmap, i0, dt)
	# drs -> distance (dx**2 + dy**2 + dz**2) between the (i+1)th and ith cluster site.
	drs.sort()
	Ns=range(1, len(drs)+1)
	#Ns.reverse()
	
	plt.figure(fnum)
	if doclf: plt.clf()
	#
	plt.loglog(drs, Ns, symb, label=label)
	plt.xlabel("$\\Delta r$")
	plt.ylabel("$N(>r)$")
	
	'''
	#dNs=[]
	Drs=[]
	for i in xrange(1, len(drs)):
		if (drs[i]-drs[i-1])==0: continue
		Drs+=[(Ns[i]-Ns[i-1])/(drs[i]-drs[i-1])]
	plt.figure(fnum+1)
	plt.loglog(Drs, Ns[1:], '-')
	'''
	return drs

def getrbsequence(Zs):
	# for now, just return rb+/rb-; we'll worry about noise, etc. later.
	rbplus=[Zs[0]]
	rbminus=[Zs[0]]
	nrbpluses=[1]
	nrbminuses=[1]
	nrbplus=1
	nrbminus=1
	#
	for z in Zs[1:]:
		rbplus+=[rbplus[-1]]
		rbminus+=[rbminus[-1]]
		#
		nrbpluses+=[nrbpluses[-1]]
		nrbminuses+=[nrbminuses[-1]]
		#print z
		if z>rbplus[-2]:
			#print "**", z
			rbplus[-1]=z
			#rbminus+=[rbminus[-1]]
			nrbplus+=1
			nrbpluses[-1]+=1
		if z<rbminus[-2]:
			#print "***", z
			rbminus[-1]=z
			#rbplus+=[rbplus[-1]]
			nrbminus+=1
			nrbminuses[-1]+=1
	#
	return [rbplus, rbminus, nrbpluses, nrbminuses]
	
def getrbs(Zs, rev=False):
	# for now, just return rb+/rb-; we'll worry about noise, etc. later.
	if rev==False:
		rbplus=Zs[0]
		rbminus=Zs[0]
		nrbplus=1
		nrbminus=1
		#
		for z in Zs[1:]:
			if z>rbplus:
				rbplus=z
				nrbplus+=1
			if z<rbminus:
				rbminus=z
				nrbminus+=1
			#
		#
		#print "nrbp: %f/%f, nrbm: %f/%f" % (nrbplus, rbplus, nrbminus, rbminus)
		return float(nrbplus)/float(nrbminus)
	if rev==True:
		rbplus=Zs[-1]
		rbminus=Zs[-1]
		nrbplus=1
		nrbminus=1
		#
		i=len(Zs)-1
		while i>=0:
			z=Zs[i]
			if z>rbplus:
				rbplus=z
				nrbplus+=1
			if z<rbminus:
				rbminus=z
				nrbminus+=1
			#
			i-=1
		#
		#print "nrbp: %f/%f, nrbm: %f/%f" % (nrbplus, rbplus, nrbminus, rbminus)
		return float(nrbminus)/float(nrbplus)

def runningBursts(zs=None, i0=0, dt=1, L=1024, geom=0, dla=1, D=1.8):
	if zs==None: zs=getNormalZs(L=L, geom=geom, dla=dla)
	zs=zs[i0:]
	cmap=ipc.getClustermap()
	drs=getDeltaRs(cmap, i0, dt)
	#
	# the basic strategy here is to dynamically define a "cluster", and the question is:
	# is this new element part of the previous cluster (which is like a branch -- and in fact,
	# explicit branching might be the better way to do this). for now, our metric will be something
	# like: if (r_corr<=n**1/D): iscluser=True
	# so the exponent will be like .5 < beta < 1.0
	# aka, if we're adding "local" sites, it's part of a cluser; we reset n when we get iscluster=False.
	#
	burstses=[]	# and w'll populate with "bursts" of [[r,z]] pairs (or we might make a class so burst.rs, burst.zs)
	#
	beta=1.0/D
	#
	thisburst=[]
	n=1
	for i in xrange(len(zs)):
		# which should be the same as len(drs)...
		if drs[i]>(n**beta):
			# "non-local"
			burstses+=[thisburst]
			thisburst=[]
			n=1
		thisburst+=[[drs[i], zs[i]]]
		n+=1
	#
	#
	return burstses

def omoribursts(burstses=None, fn0=0, nsets=5):
	# can we find something like Omori's law in these bursts? let's grab the largest clusters
	#
	if burstses==None: burstses=runningBursts()
	#
	myburstses=scipy.array(burstses).copy().tolist()
	myburstses.sort(key=lambda x:len(x))	# shortest to longest.
	myburstses.reverse()							# now longest to shortest...
	#
	for i in xrange(nsets):
		f=plt.figure(fn0+i)
		f.clf()
		#ax=gca()
		axdn=f.add_axes([.1, .05, .85, .4])
		axup=f.add_axes([.1, .5, .85, .45], sharex=axdn)
		axdn.set_xscale('linear')
		axdn.set_ylabel('cor. len $\\zeta$')
		axdn.set_yscale('linear')
		axup.set_xscale('linear')
		axup.set_yscale('linear')
		axup.set_ylabel('$z$ vals')
		
		#
		Rs=map(operator.itemgetter(0), myburstses[i])		
		Zs=map(operator.itemgetter(1), myburstses[i])
		#
		#axdn.plot(Rs, '-', label='cor. len $\\zeta$')
		axdn.bar(range(len(Rs)), Rs, 1.0, label='cor. len $\\zeta$', linewidth=0)
		axup.bar(range(len(Rs)), Zs, 1.0, label='cor. len $\\zeta$', linewidth=0)
		#axup.plot(Zs, '-', label='$z$ vals')
		X=range(1, len(Rs)+1)
		axdn.plot(X, scipy.array(X)**(1.0/1.8), '-')
		axdn.plot(X, scipy.array(X)**(1.0/2.0), '-')
	
	return myburstses
	
	
def runningburstplots(burstses=None, fn0=0, D=1.82):
	#
	if burstses==None: burstses = runningBursts()
	#
	burst_lens=scipy.array(map(len, burstses))
	burst_lens0=scipy.array(map(len, burstses))
	#logburst_sizes=map(math.log10, burst_lens)
	#
	#logburst_sizes.sort()
	burst_lens.sort()
	#
	fn=fn0
	plt.figure(fn)
	plt.clf()
	ax0=plt.gca()
	ax0.set_xscale('log')
	ax0.set_yscale('log')
	Ns=range(1, len(burst_lens)+1)
	Ns.reverse()
	#xmin=500
	#
	lfx=map(math.log10, burst_lens)
	lfy=map(math.log10, Ns)
	lfx2=[]
	lfy2=[]
	l25=math.log10(30.)
	for i in xrange(len(lfx)):
		if lfx[i]<l25: continue
		lfx2+=[lfx[i]]
		lfy2+=[lfy[i]]
	print "lens: %d, %d" % (len(lfx2), len(lfy2))
	lf1=lfp.linefit([lfx, lfy])
	plsq=lf1.doFit()[0]
	lf2=lfp.linefit([lfx2, lfy2])
	plsq2=lf2.doFit()[0]
	#plt.figure(21)
	#plt.plot(lfx, lfy, '.')
	#plt.plot(lfx, plsq[0]+scipy.array(lfx)*plsq[1], '-')
	#plt.figure(fn)
	print plsq, plsq2
	ax0.plot(burst_lens, Ns, '.-')
	X=[burst_lens[0], burst_lens[-10]]
	#a=4.0
	a=plsq[0]+.3
	b=-1.0
	#print X,  [10**(plsq[0]+plsq[1]*X[0]), 10**(plsq[0]+plsq[1]*X[1])]
	#
	#ax0.plot(X, [10**(plsq[0]+plsq[1]*X[0]), 10**(plsq[0]+plsq[1]*X[1])], 'o-', lw=2)
	ax0.plot(burst_lens, (10**plsq[0])*scipy.array(burst_lens)**plsq[1], '-', label='b=%.2f' % plsq[1])
	ax0.plot(burst_lens, (10**plsq2[0])*scipy.array(burst_lens)**plsq2[1], '-', label='b=%.2f' % plsq2[1])
	ax0.plot(burst_lens, (10**a)*scipy.array(burst_lens)**b, '--', label='b=%.2f' % b)
	#ax0.plot(X, [10**4.0, 10**4.0])
	plt.xlabel('burst length $N_{burst}$')
	plt.ylabel('Cumulative number N')
	plt.legend(loc=0)
	#
	# now plot max, min, mean z-vals (an r-vals?)
	maxZs, maxZs2, minZs, meanZs, medZs=[], [], [], [], []
	maxRs, maxRs2, minRs, meanRs, medRs=[], [], [], [], []
	meanZs2 = []
	
	geomeanZs, geomeanRs = [], []
	burst_mags=[]
	terminalZs, rbZs, rbZs2 = [], [], []
	#nup,ndown=0,0
	rbpluses, rbminuses, nrbpluses, nrbminuses=[],[],[],[]
	ts=[0]
	zs=[]
	for rw in burstses:
		# there is probably a "map" way of doing this...
		thisRs=map(operator.itemgetter(0), rw)
		thisZs=map(operator.itemgetter(1), rw)
		#
		#ts+=thisZs
		for t in thisZs:
			ts+=[ts[-1]+t]
			zs+=[t]
		#ts.pop(0)
		#
		if len(thisRs)>1:
			maxRs2 +=[max(thisRs[1:])]
			maxZs2 +=[max(thisZs[1:])]
			meanZs2 +=[scipy.mean(thisZs[1:])]
		else:
			maxRs2+=[None]
			maxZs2 +=[None]
			meanZs2 +=[None]
		maxRs  +=[max(thisRs)]
		minRs  +=[min(thisRs)]
		meanRs +=[scipy.mean(thisRs)]
		medRs  +=[scipy.median(thisRs)]
		geomeanRs +=[scipy.stats.gmean(thisRs)]
		#
		maxZs  +=[max(thisZs)]
		minZs  +=[min(thisZs)]
		meanZs +=[scipy.mean(thisZs)]
		medZs  +=[scipy.median(thisZs)]
		geomeanZs +=[scipy.stats.gmean(thisZs)]
		burst_mags += [sum(thisZs)]
		terminalZs += [thisZs[-1]]
		#
		rbZs += [getrbs(thisZs)]
		rbZs2 += [getrbs(thisZs, rev=True)]
		#if rbZs[-1]>1: nup+=1
		#if rbZs[-1]<1: ndown+=1
		#
		rbZ_sequences=getrbsequence(thisZs)
		rbpluses+=rbZ_sequences[0]
		rbminuses+=rbZ_sequences[1]
		nrbpluses+=rbZ_sequences[2]
		nrbminuses+=rbZ_sequences[3]
	#print "nup=%d, ndown=%d (%f)" % (nup, ndown, float(nup)/float(ndown))
	ts.pop(0)
	lrbZs=map(math.log10, rbZs)
	lrbZsplus=[]
	for z in lrbZs:
		if z>0:
			lrbZsplus+=[z]
		else:
			lrbZsplus+=[0.0]
	#
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	axdn=f.add_axes([.1, .05, .85, .29])
	axup=f.add_axes([.1, .33, .85, .29], sharex=axdn)
	axupup=f.add_axes([.1, .66, .85, .3], sharex=axdn)
	axupup2=axupup.twinx()
	axup.set_ylabel('z')
	axdn.set_ylabel('corlen')
	axupup.set_ylabel('burst lens')
	
	
	axdn.set_xscale('linear')
	axdn.set_yscale('linear')
	axup.set_xscale('linear')
	axup.set_yscale('linear')
	axupup.set_xscale('linear')
	axupup.set_yscale('log')
	axupup.set_title('burst time series, z,r values')
	#
	axup.plot(maxZs, 'b.-', alpha=.5, zorder=4)
	#axrb=axup.twinx()
	#axrb.plot(map(math.log10,rbZs), 'r.-', alpha=.3)
	#axrb.plot(lrbZsplus, 'r.-', alpha=.3)
	
	#axup.plot(meanZs, 'm.-', alpha=.5, zorder=5)
	axup.plot([1, len(maxZs)], [.5,.5], '-', lw=2, alpha=.7)
	axdn.plot(maxRs, 'r.-', alpha=.5)
	axdn.plot([0, len(maxRs)], [35, 35], '-')
	#axdn.plot(minRs, 'b.', alpha=.5)
	axupup.plot(burst_lens0, 'b.-', alpha=.75)
	axupup2.plot(burst_mags, 'm.-', alpha=.65, lw=2)
	#
	# expand clusters, max(Z) to full sequence length.
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	axdn=f.add_axes([.1, .05, .85, .45])
	axup=f.add_axes([.1, .55, .85, .45], sharex=axdn)
	axup2=axup.twinx()
	#axup=f.add_axes([.1, .55, .85, .45])
	#t=[]	# for now, assume interval \delta t = z
	y1=[]
	y2=[]
	y3=[]
	for i in xrange(len(burstses)):
		for j in xrange(len(burstses[i])):
			#t+=[burstses[i][1][j]]
			y1+=[maxZs[i]]
			y2+=[minZs[i]]
			y3+=[meanZs[i]]
			#
		#
	#
	#print "lens: %d, %d, %d, %d" % (len(ts), len(y1), len(y2), len(y3))
	#axup.plot(ts, y1, '-', label='maxZ', zorder=2, lw=2, alpha=.8)
	#axup.plot(ts, y2, '-', label='minZ', zorder=1, lw=1, alpha=.6)
	#axup.plot(ts, y3, '-', label='meaZ', zorder=1, lw=1, alpha=.6)
	#axup.set_xlabel('time $t$')
	axup.plot(rbpluses, 'b-', zorder=2, alpha=.95, lw=2)
	axup.plot(zs, '.', zorder=1, alpha=.5)
	#axup.plot(rbminuses,'-')
	#axup2.semilogy(scipy.array(rbpluses)/scipy.array(rbminuses), 'r-', alpha=.4)
	axup2.semilogy(scipy.array(map(float,nrbpluses))/scipy.array(map(float,nrbminuses)), 'r-', alpha=.4)
	#
	axdn.plot(y1, '-', label='maxZ', lw=1, zorder=1, alpha=.8)
	axdn.plot(zs, '.-', zorder=1, alpha=.25)
	#axdn.plot(y3, '-', label='meanZ', zorder=4, alpha=.6)
	#axdn.plot(y2, '-', label='minZ', zorder=3, alpha=.6)
	axdn.set_xlabel('natural time $n$')
	plt.legend(loc=0)
	#
	
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	plt.title('corr lens')
	ax1=plt.gca()
	ax1.set_xscale('log')
	ax1.set_yscale('log')
	#ax1.plot(burst_lens, maxRs, '.-', label='maxR')
	ax1.plot(burst_lens0, maxRs, '.', label='maxR')
	ax1.plot(burst_lens0, maxRs2, '.', label='maxR2')
	#ax1.plot(burst_lens0, minRs, '.', label='minR')
	ax1.plot(burst_lens, scipy.array(burst_lens)**(1.0/D), '-')
	ax1.set_xlabel('cluster length')
	ax1.set_ylabel('corr len $\\eta$')
	plt.legend(loc=0, numpoints=1)
	#
	#ax2=ax1.twinx()
	# mean-mean z-val:
	#meanZs2=scipy.array(meanZs2)
	#meanmean = scipy.sum((scipy.array(map(float, burst_lens0))**.5)*meanZs2)/(scipy.sum(scipy.array(map(float, burst_lens0))**.5))
	meanmean=0
	W=0
	for i in xrange(len(meanZs2)):
		#if burst_lens0[i]<100: continue
		if meanZs2[i]==None: continue
		w=float(burst_lens0[i])**.5
		meanmean+=w*meanZs2[i]
		W+=w
	meanmean/=W
	print "meanmean: %f" % meanmean
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	plt.title('Z-vals')
	ax1=plt.gca()
	ax1.set_xscale('log')
	ax1.set_yscale('linear')
	plt.plot([min(burst_lens0), max(burst_lens0)], [meanmean, meanmean], '-', lw=2, zorder=7)
	plt.plot([min(burst_lens0), max(burst_lens0)], [.25, .25], '-', lw=2, zorder=7)
	ax1.plot(burst_lens0, (scipy.array(maxZs)), '.', label='maxZ', zorder=3)
	ax1.plot(burst_lens0, (scipy.array(maxZs2)), '.', label='maxZ2', zorder=4)
	ax1.plot(burst_lens0, (scipy.array(meanZs2)), '.', label='meanZ2', alpha=.3, zorder=5)
	#ax1.plot(burst_lens0, terminalZs, '.', label='terminalZ')
	#ax1.plot(burst_lens, meanZs, '-', label='meanZ')
	#ax1.plot(burst_lens, medZs, '--', label='medZ')
	#ax1.plot(burst_lens0, minZs, '.', label='minZ')
	ax1.set_xlabel('burst length')
	ax1.set_ylabel('z')
	plt.legend(loc=0, numpoints=1)
	#
	
	'''
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	plt.title('rbZs')
	ax1=plt.gca()
	ax1.set_xscale('linear')
	#ax1.set_yscale('linear')
	#ax1.set_xscale('log')
	ax1.set_yscale('log')
	#
	#ax1.plot(burst_lens0, scipy.array(maxZs)/scipy.array(minZs), '.', label='maxZ/minZ')	
	#ax1.plot(burst_lens0, scipy.array(maxZs)/scipy.array(meanZs), '.', label='maxZ/meanZ')
	#ax1.plot(burst_lens0, terminalZs, '.')
	#ax1.plot(burst_lens0, map(math.log10,rbZs), '.')
	#ax1.plot(map(math.log10,rbZs), '.')
	ax1.plot(rbZs, '-')
	plt.legend(loc=0)
	ax1.set_xlabel('rbZs')
	#ax1.set_ylabel('z-ratio')
	ax1.set_ylabel('rb Z vals')
	'''
	
	fn+=1
	f=plt.figure(fn)
	plt.clf()
	h1=plt.hist(map(math.log10, rbZs), bins=25)
	h1=plt.hist(map(math.log10, rbZs2), bins=25, histtype='step')
	ngt, nzero, nlt=0,0,0
	for r in rbZs2:
		if r>1: ngt+=1
		if r<1: nlt+=1
		if r==1: nzero+=1
	print "rbZs2: %d, %d, %d" % (ngt, nlt, nzero)
	ngt, nzero, nlt=0,0,0
	for r in rbZs:
		if r>1: ngt+=1
		if r<1: nlt+=1
		if r==1: nzero+=1
	print "rbZs: %d, %d, %d" % (ngt, nlt, nzero)
	#ys=h1[0]
	#fn+=1
	#f=plt.figure(fn)
	#plt.clf()
	#plt.semilogy(ys, 'o-')
	
	
	return burstses
	
	
def corlenTimeSeries(zs=None, i0=0, dt=1, L=1024, geom=0, dla=1):
	if zs==None: zs=getNormalZs(L=L, geom=geom, dla=dla)
	zs=zs[i0:]
	cmap=ipc.getClustermap()
	drs=getDeltaRs(cmap, i0, dt)
	#
	plt.ion()
	f=plt.figure(0)
	plt.clf()
	#
	axdn=f.add_axes([.1, .05, .85, .4])
	axup=f.add_axes([.1, .5, .85, .4], sharex=axdn)
	#
	axdn.semilogy(drs, '-')
	axdn.set_ylabel('corlen')
	axup.plot(zs, '-')
	axup.set_ylabel('z')
	#
	rz=scipy.array(zs)*scipy.array(drs)
	'''
	plt.figure(1)
	plt.clf()
	plt.plot(rz, '-')
	'''
	#
	f4=plt.figure(4)
	plt.clf()
	ax41=f4.add_subplot(111)
	ax41.semilogy(drs, 'b-', alpha=.7, label='$\\eta$', zorder=4)
	ax41.set_xlabel('natural time $n$')
	ax41.set_ylabel('corr len')
	plt.legend(loc=0)
	#
	ax42=ax41.twinx()
	ax42.plot(zs, 'm-', alpha=.7, label='z', zorder=7)
	ax42.set_ylabel('z')
	plt.legend(loc=0)
	#
	ax43=ax41.twinx()
	ax43.plot(rz, 'c-', alpha=.7, label='$r \\times z$', zorder=6)
	ax43.set_ylabel('rz')
	#
	plt.legend(loc=0)
	#
	###############
	plt.figure(2)
	plt.clf()
	#lrzs=map(math.log10, rz)
	lrzs=map(math.log10, drs)
	plt.hist(lrzs, bins=25, log=True, normed=True)		# so with lrzs and log=True, we are getting a log-log plot of P(cor-len).
	#plt.hist(drs, bins=25, log=False)
	plt.title('P(log[$\\zeta$])')
	plt.xlabel('(log[]) corr. len')
	plt.ylabel('$p(\\zeta)$')
	
	plt.figure(3)
	plt.clf()
	ax3=plt.gca()
	ax3.set_yscale('log')
	ax3.set_xscale('log')
	#z0=.492
	z0=.5
	X=z0-scipy.array(zs)
	#X=1.0/scipy.array(zs)
	Y=scipy.array(drs)
	Xln=[min(X), max(X)]
	Yln=[min(drs), max(drs)]
	#dydx=(drs[1]-drs[0])/(X[1]-X[0])
	#print "slope: %f" % (math.log10(Yln[1]/Yln[0])/math.log10(Xln[1]/Xln[0]))
	plt.plot(X, Y, '.')
	#plt.plot(drs, X, '.')
	plt.xlabel('zs')
	plt.ylabel('corr len.')
	#
	# and a running max version of this distribution...
	plt.figure(5)
	plt.clf()
	#
	xy=[]
	for i in xrange(len(X)):
		xy+=[[X[i],Y[i]]]
	#
	xy.sort(key=lambda x: x[0])
	thisX=map(operator.itemgetter(0), xy)
	thisY=map(operator.itemgetter(1), xy)
	plt.loglog(map(operator.itemgetter(0), xy), map(operator.itemgetter(1), xy), '.', alpha=.4)
	avlen=100
	maxcorlens=[]
	for i in xrange(avlen, len(thisY)):
		maxcorlens+=[max(thisY[i-avlen:i])]
	plt.loglog(thisX[avlen:], maxcorlens, 'm.-', alpha=.5, label='avlen=%d' % avlen)
	avlen=1000
	maxcorlens=[]
	for i in xrange(avlen, len(thisY)):
		maxcorlens+=[max(thisY[i-avlen:i])]
	plt.loglog(thisX[avlen:], maxcorlens, 'c.-', alpha=.5, label='avlen=%d' % avlen)
	plt.xlabel('zs')
	plt.ylabel('corr len.')
	plt.legend(loc=0)
	#
	return zs

def getDeltaRs(cmap, i0, dt=1):
	# distance between subsequent sites.
	# assume cmap is like [ [i, x, y, z, val] ]
	X=map(operator.itemgetter(1), cmap)
	Y=map(operator.itemgetter(2), cmap)
	Z=map(operator.itemgetter(3), cmap)
	#
	drs=[]
	print "set range: ", (i0+dt), ", ", len(cmap)
	for i in xrange((i0+dt), len(cmap)):
		x2, y2, z2 = X[i], Y[i], Z[i]
		x1, y1, z1 = X[i-dt], Y[i-dt], Z[i-dt]
		#
		dr = math.sqrt((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)
		drs+=[dr]
		#
	#
	return drs
		

def occProbL(nits=100, lw0=4, lwmax=10, geom=0, dla=1):
	# occupation probability (aka, N/(L*L)
	# we want <cluster size>(L) for 2D lattice(s).
	# to be a little bit schmancy, let's pick log(L) randomly and just take some number of samples. we can use arbitrary L, but
	# it seems right to stick to powers of 2.
	#nits=1000
	R1=random.Random()
	L,N,Nb = [],[],[]
	Nave, Nbave = [],[]
	
	print "OccProb for: nits=%d, log(width) = (%d - %d), geom=%d, dla=%d" % (nits, lw0, lwmax, geom, dla)
	#
	width=2**lw0
	
	percobj=ipc.ipercObj(int(width), geom, dla)
	rseed=ipc.getrandSeed()
	#
	for j in xrange(lw0, lwmax):
		width=2**j
		ipc.init(width, geom, rseed)
		N0, Nb0=[],[]
		for i in xrange(nits):
			#width=2**R1.randint(lw0, lwmax)
			#csize=ipc.ipercObj(int(width), geom, dla)
			#ipc.init(width, geom, rseed)
			#print "initted..."
			ipc.initializeLattice()
			ipc.initializeCluster()
			csize=ipc.makePercolatingCluster()
			csize2=ipc.getClustsize()
			csize3=len(ipc.getClusterBonds())
			#
			#ipc.init(width, geom, rseed)
			#
			rseed+=1
			print i, ": ", csize, csize2, csize3, "width=", width
			#
			L+=[width]
			N0+=[float(csize)/float(width*width)]
			Nb0+=[float(csize3)/float(width*width*ipc.getnNN())]
			N+=[N0[-1]]
			Nb+=[Nb0[-1]]
		#N+=[N0]
		#Nb+=[Nb0]
		Nave+=[scipy.mean(N0)]
		Nbave+=[scipy.mean(Nb0)]
	#
	print len(L), len(N)
	plt.ion()
	plt.figure(0)
	plt.clf()
	plt.loglog(L,N, 'bo', alpha=.35)
	plt.loglog(L,Nb, 'go', alpha=.35)
	#
	#plt.figure(1)
	#plt.clf()
	print set(L)
	print Nave
	L1=list(set(L))
	L1.sort()
	plt.loglog(L1, Nave, 'bo-')
	plt.loglog(L1, Nbave, 'go-')
	
	return [L, N, Nb]
	
def plotClust(xy=[], fnum=0):
	# assume from getClusterMap()
	X=map(operator.itemgetter(1), xy)
	Y=map(operator.itemgetter(2), xy)
	#
	plt.ion()
	plt.figure(fnum)
	plt.plot(X,Y, '.')
			
def getXYfromClust(Zs, width):
	X1, Y1 = [], []
	#
	for rw in Zs:
		X1+=[rw%width]
		Y1+=[rw/width]
	#
	return [X1, Y1]
##
def getPerc2(nits=1, w=64, geom=3, dla=1):
	csize=ipc.ipercObj(w, geom, dla)
	#
	#csize=ipc.getClustsize();
	#ipc.initializeLattice()
	ipc.initializeCluster()
	thissize=ipc.makePercolatingCluster()
	print thissize
	csize=ipc.getClustsize()
	#ks+=[csize]
	ipc.randSeedPlus()
	#
	plt.ion()
	fig1=plt.figure(1)
	plt.clf()
	ax1 = fig1.gca(projection='3d')
	#
	xyz = ipc.getClustermap()
	X=map(operator.itemgetter(1), xyz)
	Y=map(operator.itemgetter(2), xyz)
	Z=map(operator.itemgetter(3), xyz)
	#
	ax1.plot(X,Y,Z, '.')
	#
	'''
	fig2=plt.figure(2)
	plt.clf()
	plt.ion()
	ax2 = fig2.gca(projection='3d')
	#
	Zs = ipc.getClusterSites()
	X2,Y2,Z2=[],[],[]
	#
	for rw in Zs:
		X2+=[rw%w]
		Y2+=[(rw/w)%w]
		Z2+=[rw/(w*w)]
		
		#print X2[-1], Y2[-1], Z2[-1]
	#
	ax2.plot(X2,Y2,Z2, '.')
	'''
	
	#return [X,Y,Z]
	return xyz
	
	
	
# Development and diagnostic:
def getPerc1(width=256, geom=0, dla=1, i0=None):
	csize=ipc.ipercObj(width, geom, dla)
	if i0==None:
		#i0=(width*width + width)/2
		#if geom==3: i0+=(width*width*width)/2
		q=ipc.initializeCluster()
	else:
		q=ipc.initializeCluster(i0)
	#
	
	thissize=ipc.makePercolatingCluster()
	rs1=ipc.getrandSeed()
	Z=ipc.getClusterSites()
	Z2=ipc.getLatticeSites()
	#	
	X1, Y1, X2, Y2 = [], [], [], []
	#
	for rw in Z:
		X1+=[rw%width]
		Y1+=[rw/width]
	#
	for i in xrange(len(Z2)):
		if Z2[i]!=1: continue
		#print (i-1), (i-1)%256, (i-1)/256
		X2+=[(i-0)%width]
		Y2+=[(i-0)/width]
	#
	print csize, len(X1), len(X2)
	print Z[0:5]
	print Z2[0:5]
	#
	plt.ion()
	plt.figure(0)
	plt.clf()
	plt.plot(X1, Y1, '.')

	plt.figure(1)
	plt.clf()
	plt.plot(X2, Y2, '.')
	#
	#####
	ipc.randSeedPlus()
	csize=ipc.ipercObj(width, geom, dla, rs1+1)
	thissize=ipc.makePercolatingCluster()
	Z=ipc.getClusterSites()
	Z2=ipc.getLatticeSites()
	
	X1, Y1, X2, Y2 = [], [], [], []
	#
	for rw in Z:
		X1+=[rw%width]
		Y1+=[rw/width]
	#
	for i in xrange(len(Z2)):
		if Z2[i]!=1: continue
		#print (i-1), (i-1)%256, (i-1)/256
		X2+=[(i)%width]
		Y2+=[(i)/width]
	#
	plt.ion()
	plt.figure(2)
	plt.clf()
	plt.plot(X1, Y1, '.')

	plt.figure(3)
	plt.clf()
	plt.plot(X2, Y2, '.')
	#
	#return [X1, Y1, X2, Y2]
	return Z2

def getPerc3(w=16, geom=0, dla=1):
	# test percolation model:
	#1) get a wxw lattice
	#2) get sites and edges
	#3) write both to file
	#4) makePercolatingCluster
	#5) compare.
	# note: we can also go step by step.
	a=ipc.init(w, geom)
	nnn=ipc.getnNN()
	print "nnn? " , nnn
	#if geom==0 or geom==1: nnn=4
	#if geom==2 or geom==3: nnn=6
	#print "nnn: %d" % nnn
	#
	
	'''
	sites=ipc.getLatticeSites()
	edges=ipc.getLatticeBonds()
	fout=open('data/sites.dat', 'w')
	for i in xrange(w):
		for j in xrange(w):
			fout.write('%d\t' % sites[w*i+j])
		fout.write('\n')
	fout.close()
	#
	fout=open('data/edges.dat', 'w')
	for i in xrange(w):
		for j in xrange(w):
			fout.write("%d:" % (w*i + j))
			for k in xrange(nnn):
				fout.write('%d,' % (edges[w*nnn*i + nnn*j + k][2]))
			fout.write('\t')
		fout.write('\n')
	fout.close()
	fout=open('data/edges2.dat', 'w')
	for i in xrange(len(edges)):
		fout.write('%d\t%d\t%d\n' % (edges[i][0], edges[i][1], edges[i][2]))
	fout.close()
	#
	'''
	pyclusts=[]
	pyclusts+=[ipc.initializeCluster()]
	#print pyclusts
	site0=ipc.getClusterSites()[0]
	#
	csize=ipc.makePercolatingCluster(dla)
	clustE=ipc.getClusterBonds()
	
	plotClustBonds(clustE, fnum=0, w=w, geom=geom, dla=dla)
	# now, let's look at strands:
	edges=ipc.getClusterBonds()
	#print edges[0:10]
	strands = getStrands(edges)
	#print strands[0:3]
	
	#
	#print "bonds: %d, strands: %d" % (len(edges), len(strands))
	#
	plotStrands(strands)
	#
	return None

def plotTipSites(fnum=0, doclf=False):
	tips=ipc.getTips()
	ends = map(operator.itemgetter(2), tips)
	starts = map(operator.itemgetter(1), tips)
	W=ipc.getwidth()
	#geom=ipc.getgeom()	// for now assume square...
	#
	Xends=[]
	Yends=[]
	Xstarts=[]
	Ystarts=[]
	for i in xrange(len(tips)):
		Xends+=[tips[i][2]%W]
		Xstarts+=[tips[i][1]%W]
		Yends+=[tips[i][2]/W]
		Ystarts+=[tips[i][1]/W]
	#
	f1=plt.figure(fnum)
	if doclf==True: plt.clf()
	#
	plt.plot(Xends, Yends, 'bo')
	plt.plot(Xstarts, Ystarts, 'gd')
	#
	return [Xends, Yends, Xstarts, Ystarts]

def clusterMovie(clustE=None, fnum=0, w=None, geom=None, dla=1, doclf=True, clr0=0, outdir='movies/movie1', blockout=True):
	if geom==None: geom=ipc.getgeom()
	if w==None: w=ipc.getwidth()
	if w==None:
		if geom in (0,1,2): w=math.sqrt(len(clustE))
		if w==3: w=(len(clustE))**(1./3.)
	if clustE==None:
		a=ipc.init(w, geom)
		ipc.setdlamode(dla)
		b=ipc.initializeCluster()
		c=ipc.makePercolatingCluster()
		clustE=ipc.getClusterBonds()
	#
	
	nnn=ipc.getnNN()
	plt.figure(fnum)
	plt.clf()
	#strlen=int(math.log10(len(clustE)) + 1)
	strlen=7
	framestrnum=''
	for i in xrange(strlen):
		framestrnum=framestrnum+'0'
		
	for i in xrange(1,len(clustE)):
		if blockout==True:
			plt.plot([0, w], [0,w], ',', alpha=.01)
		#plotClustBonds(clustE[0:i], fnum=fnum, w=w, geom=geom, dla=dla, doclf=True, clr0=clr0)
		plotClustBonds(clustE[i-1:i], fnum=fnum, w=w, geom=geom, dla=dla, doclf=False, clr0=clr0)
		numstr=(framestrnum + str(i))[-strlen:]
		plt.savefig('%s/img%s.png' % (outdir, numstr))
	#
	avconvstr = 'avconv -i %s/img%07.png -r 50 %s/iperc-%d-%d.mp4' % (outdir, outdir, w, geom)
	#avconvstr='avconv -i %scontour/%s-%05d.png -r 50 %s/%sconts.mp4' % (moviedir, nameroot, moviedir, nameroor)
	os.system(avconvstr)
	
def plotClustBonds(clustE=None, fnum=0, w=None, geom=None, dla=1, doclf=True, clr0=0):
	if clustE==None: clustE=ipc.getClusterBonds()
	if geom==None: geom=ipc.getgeom()
	if w==None: w=ipc.getwidth()
	if w==None:
		if geom in (0,1,2): w=math.sqrt(len(clustE))
		if w==3: w=(len(clustE))**(1./3.)
	#
	colors = {0:'b', 1:'g', 2:'r', 3:'c', 4:'m', 5:'y', 6:'k'}
	site0=ipc.getClusterSites()[0]
	nnn=ipc.getnNN()
	fig0=plt.figure(fnum+0)
	if doclf: plt.clf()
	clr0=clr0%len(colors)
	clr=clr0
	#
	if geom in (0, 1, 2):
		x0=site0%w
		y0=site0/w
		plt.plot([x0], [y0], '*', ms=15)
		#	
		#Xf,Yf, Xt, Yt=[], [], [], []
		#for rw in clustE:
		#
		#ncolors=0
		for i in xrange(len(clustE)):
			rw = clustE[i]
			lbl=None
			if i==0: lbl='(0)'
			if i>0 and (rw[0]!=clustE[i-1][1]):
				# new branch.
				clr0+=1
				clr=(clr0)%(len(colors))
				#
				if (clr0)/len(colors) == 0:
					lbl = '(%d)' % ((clr))
				
			#Xf+=[rw[0]%w]
			#Xt+=[rw[1]%w]
			#Yf+=[rw[0]/w]
			#Yt+=[rw[1]/w]
			dx0=0.
			dx1=0.
			y0=rw[0]/w
			y1=rw[1]/w
			if geom==2:
				dx0=+.5*(y0)%2.
				dx1=+.5*(y1)%2.
			#
			plt.plot([rw[0]%w + dx0, rw[1]%w + dx1], [rw[0]/w, rw[1]/w], '.-', color=colors[(clr)], lw=1.0 + 1.0*float(rw[2]/(w*w)), alpha=.25+.75*float(rw[2]/(w*w)), label=lbl)
	#
	if geom==3:
		# cubic, 3D geometry...
		x0=site0%w
		y0=(site0/w)%w
		z0=site0/(w*w)%w		# note: x_{j = 0,1,2} = (i/(w**j))%w
		ax=plt.gca()
		ax1 = fig0.gca(projection='3d')
		#
		# initialization point:
		plt.plot([x0], [y0], [z0], '*', ms=15)
		#	
		#Xf,Yf, Xt, Yt=[], [], [], []
		for i in xrange(len(clustE)):
			rw = clustE[i]
			lbl=None
			if i==0: lbl='(0)'
			if i>0 and (rw[0]!=clustE[i-1][1]):
				# new branch.
				clr0+=1
				clr=(clr0)%(len(colors))
				#
				if (clr0)/len(colors) == 0:
					lbl = '(%d)' % ((clr))
			#
			x0=rw[0]%w
			x1=rw[1]%w
			#
			y0=(rw[0]/w)%w
			y1=(rw[1]/w)%w
			#
			z0=rw[0]/(w*w)
			z1=rw[1]/(w*w)			
			#
			plt.xlabel('x')
			plt.ylabel('y')
			#ax.set_zlabel('z')
			#plt.plot([x0, x1], [y0, y1], [z0, z1], '--', lw=1.0 + 1.0*float(rw[2]/(w*w)), alpha=.25+.75*float(rw[2]/(w*w)))
			
			plt.plot([x0, x1], [y0, y1], [z0, z1], '.-', color=colors[clr], label=lbl)
	#plt.legend(('(1)','(2)', '(3)', '(4)', '(5)', '(6)', '(7)'), loc='best')
	plt.legend(loc='best', numpoints=2)
	#
	return None

def xyzsfromindex(rws, w=None, geom=None):
	if w==None: w=ipc.getwidth()
	if geom==None: geom=ipc.getgeom()
	#
	X0, X1, Y0, Y1, Z0, Z1 = [], [], [], [], [], []
	#
	for rw in rws:
		A=xyzfromindex(rw, w, geom)
		X0+=[A[0]]
		X1+=[A[1]]
		Y0+=[A[2]]
		Y1+=[A[3]]
		Z0+=[A[4]]
		Z1+=[A[5]]
	#	
	return [X0, X1, Y0, Y1, Z0, Z1]
	
def xyzfromindex(rw, w=None, geom=None):
	if w==None: w=ipc.getwidth()
	if geom==None: geom=ipc.getgeom()
	#
	# rw sill be and edge row: (from, to , val)
	if geom in (0,1,2):
		dx0=0.
		dx1=0.
		y0=rw[0]/w
		y1=rw[1]/w
		if geom==2:
			dx0=+.5*(y0)%2.
			dx1=+.5*(y1)%2.
		x0=rw[0]%w + dx0
		x1=rw[1]%w + dx1
		z0=0.
		z1=0.
	if geom == 3:
		x0=rw[0]%w
		x1=rw[1]%w
		#
		y0=(rw[0]/w)%w
		y1=(rw[1]/w)%w
		#
		z0=rw[0]/(w*w)
		z1=rw[1]/(w*w)
	#
	#return [[x0, y0, z0], [x1, y1, z1]]
	return [x0, x1, y0, y1, z0, z1]
#
def segstrandfmplot(clustE=None, fnum=0):
	if clustE==None: clustE=ipc.getClusterBonds()
	segs=getsegs(clustE)
	strands=getStrands(clustE)
	#
	seg_lens=map(len, segs)
	strand_lens=map(len, strands)
	Nsegs=range(1, len(segs)+1)
	Nstrands=range(1, len(strands)+1)
	#
	seg_lens.sort()
	strand_lens.sort()
	Nsegs.reverse()
	Nstrands.reverse()
	#
	f=plt.figure(fnum)
	plt.clf()
	ax=plt.gca()
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.plot(seg_lens, Nsegs, '-')
	ax.plot(strand_lens, Nstrands, '-')
	#
	return clustE
	

def getsegs(clustE=None):
	# 
	#
	#clustE -> [ [from, to, val], []...]
	if clustE==None: clustE=ipc.getClusterBonds()
	segs=[ [clustE[0]] ]
	for i in xrange(1, len(clustE)):
		thisbond=clustE[i]
		#
		if (thisbond[0]!=segs[-1][-1][1]):
			segs+=[[]]
		segs[-1]+=[thisbond]
	return segs
	
def getStrands(clustE=None):
	# note: this is geometry independend (until we want to plot them)
	# start at the seed-site; follow it until start != finish; initiate a new list.
	# check existing strands to see if we're picking up off the end of one, otherwise start a new list.
	# returns struct is: list-o-strands-o-edges:
	# [ [strand [edges]] ]
	# clustE is the time-ordered list of bonds.
	#
	#clustE -> [ [from, to, val], []...]
	if clustE==None: clustE=ipc.getClusterBonds()
	strands=[ [clustE[0]] ]
	strandindex=0
	#
	for i in xrange(1,len(clustE)):
		thisbond=clustE[i]
		#
		#print "this: %d, prev: %d" % (thisbond[0], strands[strandindex][-1][1])
		#print "bonds: ", thisbond, strands[strandindex][-1]
		if (thisbond[0]!=strands[strandindex][-1][1]):
			# this bond's "from" != previous bond's "to"
			# look at all previous clusters. do their end-sites match this bond?
			# make new strand and set index, in the event we don't find a new strand.
			strandindex=len(strands)
			strands+=[[]]
			for j in xrange(len(strands)-1):
				if strands[j][-1][1]==thisbond[0]:
					# the new edge matches this strand.
					strandindex=j
					strands.pop(-1)	# pop holding list we just created.
					break
				#
			#
		strands[strandindex]+=[thisbond]
	#
	return strands

def plotStrands(strands=None, fnum=1, linestring='.-'):
	# use: xyzfromindex(rw, w, geom):
	if strands==None: strands=getStrands()
	plt.figure(fnum)
	plt.clf()
	w=ipc.getwidth()
	#
	lens=[]
	#
	clrindex=0
	for strand in strands:
		#xyzs = xyzsfromindex(strand, w=None, geom=None)
		'''
		X=[]
		Y=[]
		for edge in strand:
			X+=[edge[0]%w, edge[1]%w]
			Y+=[edge[0]/w, edge[1]/w]
		#plt.plot(X,Y, linestring)
		#lens+=[len(X)]
		'''	
		#
		#plt.plot([xyzs[0][0]] + xyzs[1], [xyzs[2][0]] + xyzs[3], '.-')
		#lens+=[len(xyzs[0])]
		lens+=[len(strand)]
		# plotClustBonds(clustE, fnum=0, w=w, geom=geom, dla=dla)
		z=plotClustBonds(strand, fnum=fnum, doclf=False, clr0=clrindex)
		clrindex+=1
		
	#
		# assume normal start:
	#i=(w*w + w)/2
	#plt.plot([i%w], [i/w], 'b*', ms=15)

	lens.sort()
	lens.reverse()
	plt.figure(fnum+1)
	plt.clf()
	plt.loglog(lens, range(1, len(lens)+1))
	
	#



def plotipercFile(fname='gridout2.dat'):
	fin = open(fname, 'r');
	#
	X, Y, Z, I = [], [], [], []
	Xnull, Ynull=[],[]
	XY=[]
	datalen=0
	dx=0
	dy=0	# offsets for different lattice types
	intltype=0
	strltype='square'
	#rownum=0
	#
	for rw in fin:
		
		if "#lattype\t" in rw:
			rws=rw.split()
			print rws
			intltype = int(rws[1])
			strltype = rws[2]
			width = int(rws[3])
			continue
		#
		if rw[0]=='#': continue
		
		rws=map(int, rw.split())
		if intltype == 2:
			dx = .5*((datalen/width)%2)
			#dx=0
		datalen+=1
		if rws[3]!=1: 
			Xnull+=[rws[1]+dx]
			#Ynull+=[(math.sqrt(3)/2.0)*(rws[2] + dy)]		# something like this for tri-lattice.
			Ynull+=[rws[2] + dy]
			continue
		#	
		#
		I+=[rws[0]]
		#X+=[rws[1] + dx]
		X+=[rws[0]%256]
		#Y+=[(math.sqrt(3)/2.0)*(rws[2] + dy)]
		#Y+=[(rws[2] + dy)]
		Y+=[rws[0]/256]
		Z+=[rws[3]]
		#
		XY+=[[rws[1]+dx, rws[2]+dy]]
		
		#
	fin.close()
	#
	#datalen=len(X)
	width=math.sqrt(datalen)
	
	#
	plt.figure(0)
	plt.clf()
	plt.ion()
	#
	plt.plot([0, width+1, width+1, 0, 0], [0, 0, width+1, width+1, 0], 'k-', zorder=0)
	plt.plot(X,Y, 'b.', ms=5)
	plt.plot(Xnull, Ynull, 'g.', ms=3, alpha=.6)
	#
	plt.figure(1)
	plt.clf()
	plt.plot([0, width+1, width+1, 0, 0], [0, 0, width+1, width+1, 0], 'k-', zorder=0)
	for xy in XY:
		x = xy[0]
		y = xy[1]
		#
		X0=[x, x+1, x+1, x, x]
		Y0=[y, y, y+1, y+1, y]
		#
		plt.plot(X0,Y0, 'b-')
		plt.fill(X0,Y0, 'b', alpha=.45)
	#
	'''
	plt.figure(2)
	plt.clf()
	plt.plot([0, width+1, width+1, 0, 0], [0, 0, width+1, width+1, 0], 'k-', zorder=0)
	for xy in XY:
		x = xy[0]
		y = xy[1]
		#
		plt.plot([x], [y], 'g.', zorder=3)
		X0=[x-1, x+1]
		Y0=[y, y]
		plt.plot(X0,Y0, 'b.', zorder=2)
		#
		X0=[x, x]
		Y0=[y-1, y+1]
		plt.plot(X0,Y0, 'b.',  zorder=2)
	'''
	#
	for a in [X,Y,Z,I, XY, X0, Y0]:
		a=None

#
def intervalDistAnalysis(catname='cats/socaliperc.cat', lon=[-122.5, -114.0], lat=[30.5, 37.25], minMag=2.5, dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None], Nmax=999999, beta1=.5, beta2=1.0, refreshcat=True, dataset=None):
	# catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.datetime(2005,01,01, tzinfo=tzutc), None], Nmax=999999, fout='cats/mycat.cat')
	#tohoku:
	# lons=[135., 146.], lats=[30., 41.5]
	# approx for norcal:
	# ipp.intervalDistAnalysis(catname='cats/norcaliperc.cat', lon=[-128., -117.], lat=[37., 43.], minMag=3.0, refreshcat=True)
	if dataset=='tohoku':
		catname='cats/tohokutmp.cat'
		lon=[135., 150.]
		lat=[30., 41.5]
		minMag=4.5
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,03,01, tzinfo=pytz.timezone('UTC'))]
		dates0=[dtm.datetime(2011,3,11, 12, tzinfo=pytz.timezone('UTC')), None]
	if dataset=='socal':
		catname='cats/socaliperc.cat'
		lon=[-122.5, -114.0]
		lat=[30.5, 37.25]
		minMag=2.5
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
		dates0=[dtm.datetime(2009,01,01, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,01,01, tzinfo=pytz.timezone('UTC'))]
	if dataset=='norcal':
		catname='cats/norcaliperc.cat'
		lon=[-128., -117.]
		lat=[37., 43.]
		minMag=3.0
		dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
	#
	if refreshcat:
		clist = atp.catfromANSS(lon=lon, lat=lat, minMag=minMag, dates0=dates0, Nmax=Nmax, fout=catname)
		clist.sort(key=lambda x:x[0])
		c1=eqp.eqcatalog(clist)
	else:
		c1=eqp.eqcatalog()
		c1.loadCatFromFile(catname)
	'''
	c1=eqp.eqcatalog()
	c1.loadCatFromFile('cats/pfshock-20120109.cat')
	'''
	#
	interval_windowlen=100
	fdts=map(mpd.date2num, map(operator.itemgetter(0), c1.getcat(0)))
	fdts=scipy.array(fdts)
	intervals = fdts[1:]-fdts[0:-1]
	intervals2= fdts[interval_windowlen:]-fdts[0:-interval_windowlen]
	mags=map(operator.itemgetter(3), c1.getcat(0)[:-1])
	overNprimes=10.0**(minMag - scipy.array(mags))
	intprimes=intervals*overNprimes
	#print "**", intprimes[0:10]
	#print "**", intervals[0:10]
	intervalsesprime=intprimes[1:]-intprimes[:-1]
	#rates = 1./scipy.array(intervals)
	'''
	rates=[]
	for dt in intervals:
		#print dt
		if dt!=0.0:
			rates+=[1.0/dt]
			#print "rate/interval: %f, %f" % (rates[-1], dt)
	#
	rates=scipy.array(rates)
	rateses = rates[1:] - rates[:-1]
	print "ratelens: %d, %d" % (len(rateses), len(rates))
	'''
	
	intervalses = intervals[1:] - intervals[:-1]
	intervalses2 = []
	intsGT=[]
	intsLT=[]
	for i in xrange(1, len(intervals)):
		dt2=intervals[i]
		dt1=intervals[i-1]
		#print dt
		if dt1!=0.0 and dt2!=0:
			intervalses2+=[dt2/dt1]
			#print "rate/interval: %f, %f" % (rates[-1], dt)
		#
		ddt=dt2-dt1
		if ddt>0: intsGT+=[ddt]
		if ddt<0: intsLT+=[-ddt]
	#	
	#
	
	plt.figure(0)
	plt.clf()
	plt.ion()
	ax=plt.gca()
	ax.set_yscale('log')
	#ax.plot(fdts[interval_windowlen:], intervals2, '-')
	reversed_intervals = scipy.array(intervals2).copy().tolist()
	reversed_intervals.reverse()
	#ax.plot(fdts[interval_windowlen:], reversed_intervals, '-')
	ax.plot(intervals2, '.')
	ax.plot(reversed_intervals, '-')
	#ax2=ax.twinx()
	#ax2.plot(fdts[interval_windowlen:], 1.0/scipy.array(intervals2), 'r-')
	#plt.plot(fdts[1:], intprimes, '.')
	#c1.plotCatMap(fignum=0)
	#
	#return None
	
	f1=plt.figure(1)
	plt.clf()
	ax1=f1.gca()
	ax2=ax1.twinx()
	plt.title('interval distribution, $t_i - t_{i-1}$')
	h1=ax1.hist(intervalses, bins=200, normed=True, log=True, label='intervals')
	#lN1 = map(math.log10, h1[0])
	h3=ax1.hist(intervalsesprime, bins=200, normed=True, log=True, histtype='step', label='intervals/n(m)')
	#h3=ax1.hist(intervalsesprime[1:]-intervalsesprime[:-1], bins=200, normed=True, log=True, histtype='step')
	plt.legend(loc=0)
	#
	#hdiff=scipy.array(h1[0])/scipy.array(h3[0])
	#ax2.plot(h3[1][1:], hdiff, 'r-', lw=2)
	#
	plt.figure(2)
	plt.clf()
	#print intervalses2[0:10]
	lints=map(math.log10, intervalses2)
	#h2=plt.hist(rateses, bins=100, normed=True, log=True)
	plt.title('log(log[]) distribution')
	h2=plt.hist(lints, bins=200, normed=True, log=True)
	#
	#plt.figure(3)
	#plt.clf()
	#lints2=map(math.log10, intervalsesprime)
	#
	# so maybe z != \delta t, which is not too shocking since \deltat t>0 and 0<z<1.0(.5).
	# what if we start with a Poisson distribution:
	# p(dt) = 1-exp(-dt/<dt>), and for random numbers et al, p->u, and u=1-u, AND u -> z. so,
	# z = exp(-dt/<dt>), or maybe we use a PL distribution
	# we can either get
	#
	f=	plt.figure(3)
	plt.clf()
	ax=plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('linear')
	ax.hist(intsGT, bins=100, histtype='step')
	ax.hist(intsLT, bins=100, histtype='step')
	plt.title('$z$ distribution(s)')
	#
	zs0 = getZsfromInts(intervals, zfact=1.0, beta=beta1)
	#beta2=.5*betaweib
	zs1 = getZsfromInts(intervals, zfact=1.0, beta=(beta2))
	#zs1 = getZsfromInts(intprimes, zfact=1.0, beta=betaweib)
	zstrunc1=[]
	zthresh=.5
	for z in zs0:
		if z>zthresh: zstrunc1+=[z-zthresh]
	#
	f=plt.figure(4)
	plt.clf()
	ax=plt.gca()
	plt.title('Hist of zs[i]-zs[i-1]')
	plt.ylabel('$p(\\Delta z_1)$')
	plt.xlabel('$\\Delta z_1$')
	#ax.set_yscale('log')
	#ax.set_xscale('linear')
	dzs=zs0[1:]-zs0[:-1]
	dzs1=zs1[1:]-zs1[:-1]
	h0=ax.hist(dzs, bins=100, normed=True, log=True)
	h1=ax.hist(dzs1, bins=100, normed=True, log=True, histtype='step')
	#
	GT,LT=[],[]
	#for d in dzs:
	#	if d>0: GT+=[d]
	#	if d<0: LT+=[-d]
	f=plt.figure(5)
	plt.clf()
	plt.title('Split Hist of zs[i]-zs[i-1]')
	#plt.hist(GT, bins=50, normed=True, log=True, histtype='step')
	#plt.hist(LT, bins=50, normed=True, log=True, histtype='step')
	Hlt=(h0[0][0:len(h0[0])/2]).tolist()
	Hlt.reverse()
	ax=plt.gca()
	ax.set_xscale('linear')
	ax.set_yscale('log')
	plt.plot(Hlt, '-', label='lts')
	plt.plot(h0[0][len(h0[0])/2:], '-', label='gts')
	plt.ylabel('$p(\\Delta z_1)$')
	plt.xlabel('$\\Delta z_1$')
	plt.legend(loc=0)
	#
	f=plt.figure(6)
	plt.clf()
	ax=plt.gca()
	ax.hist(zs0, bins=200, normed=True, log=True, label='$\\beta=%.2f$' % beta1, histtype='step', zorder=5)
	ax.hist(zs1, bins=200, normed=True, log=True, label='$\\beta=%.2f$' % beta2, histtype='step', zorder=4)
	ax.hist(.5*(scipy.array(zs0)+scipy.array(zs1)), bins=200, normed=True, log=True, label='$(zs0+zs1)/2$', zorder=3)
	ax.hist(zstrunc1, bins=200, normed=True, log=True, label='h1-truncated', histtype='step', zorder=6)
	plt.legend(loc=0)
	plt.title('$z$ distribution, Weibull transformation with $\\beta _1 = %.2f$, $\\beta _2 = %.2f$' % (beta1, beta2))
	#
	f=plt.figure(7)
	plt.clf()
	ax=plt.gca()
	nbins=50
	#ax.set_xscale('linear')
	#ax.set_yscale('log')
	h1=ax.hist(scipy.array(zstrunc1[1:])-scipy.array(zstrunc1[:-1]),bins=nbins, normed=True, log=False, label='z1-truncated', histtype='bar', zorder=3)
	ax.plot([0.,0.], [min(h1[0]), 1.1*max(h1[0])], 'r-', lw=2, alpha=.7, zorder=5)
	lt=h1[0][0:nbins/2]			# "less than's"
	ltx=-1.0*scipy.array(h1[1][0:nbins/2])
	gt=h1[0][nbins/2:]
	gtx=h1[1][nbins/2+1:]
	plt.plot(ltx, lt, '-', zorder=6, label='lt')
	plt.plot(gtx, gt, '-', zorder=7, label='gt')
	plt.title('$z_i$ truncated, markov dist.')
	plt.legend(loc=0)
	#
	f=plt.figure(8)
	plt.clf()
	ax=plt.gca()
	lt=lt.tolist()
	lt.reverse()
	diffs=(scipy.array(gt)-scipy.array(lt))
	ratios=(scipy.array(gt)/scipy.array(lt))
	#ratios=(scipy.array(gt)/scipy.array(lt))**(-scipy.array(gt))
	#ratioslt=(scipy.array(gt)/scipy.array(lt))/scipy.array(map(math.log10, lt))
	#ratios=(scipy.array(map(math.log10, gt)) - scipy.array(map(math.log10, lt)))/scipy.array(map(math.log10, gt))
	diffsum=[0.]
	#diffgeom=[1.0]
	gtltGeom=[]		# running geom-mean.
	glltmeans=[]
	i=0
	for d in diffs:
		diffsum  +=[diffsum[-1]+d]
		#diffgeom +=[diffgeom[-1]*d]	# this gets really small, really fast, so it's probably not useful.
		if i>=0:
			gtltGeom+=[scipy.stats.gmean(ratios[0:i+1])]
			glltmeans+=[scipy.mean(ratios[0:i+1])]
		else:
			gtltGeom+=[1.0]
		i+=1
	#
	#ax.plot(scipy.zeros(len(diffs)), '-')
	ax.plot([min(gtx), max(gtx)], [0., 0.], 'k-')
	ax.plot(gtx,diffs, 'b-', label='(p(gt)-p(lg))')
	ax.plot(gtx,diffsum[1:], 'r-', label='$\\sum {(p[gt]-p[lt])}$')
	plt.legend(loc=0)
	#
	#ax2=ax.twinx()	
	#ax2.plot([min(gtx), max(gtx)], [1., 1.], 'k-')
	#ax2.plot(gtx, gtltGeom, 'c-')
	#ax2.plot(gtx, glltmeans, 'y-')
	#ax2.plot(gtx, ratios, 'm-')
	#ax2.plot([min(gtx), max(gtx)], scipy.ones(2)*scipy.mean(ratios), 'k--')
	#
	#ax2.plot(gtx, diffgeom[1:], 'r.--')
	#ax2.plot(gtx, ratioslt, 'y-')
	#ax2.plot(gtx, 1.0/scipy.array(map(math.log10, gt)), 'c--')
	#ax2.plot(gtx, 1.0/scipy.array(map(math.log10, lt)), 'r--')
	#ax2.plot([min(gtx), max(gtx)], [1., 1.], 'k-')
	#ax2.plot([min(gtx), max(gtx)], [0., 0.], 'k-')
	ax1.legend(loc=0)
	ax2.legend(loc=0)
	#
	print "gtX: ", gtx[0:10]
	print "data: ", gt[0:10], lt[0:10], diffs[0:10]
	#
	return c1
	
	
def getIntsfromZs(Zs, beta=.5, tbar=1.0):	
		deltat = (-tbar*scipy.array(map(math.log10, Zs)))**(1.0/beta)
		# noting that 0<z<1, and i guess tbar (average t or otherwise the exponential denominator) is arbitrary.
		#
		return deltat
	
def getZsfromInts(ints, zfact=1.0, beta=.5):
	# extraxt z-values from intervals. for now, assume poisson distribution.
	# intervals are all positive, so use a regular mean.
	tbar=scipy.mean(ints)
	# poisson:
	#zs=1.0-scipy.array(map(math.exp, -scipy.array(ints)/tbar))
	#
	# what about a Weibull distribution?
	#beta=0.5	# implying that prob. of failure increases in time.
	#zs=1.0-scipy.array(map(math.exp, -(scipy.array(ints)/tbar)**beta))
	zs=fweibull(x=scipy.array(ints), lamb=tbar, beta=beta)
	#
	return zs*zfact
	
def fweibull(x, lamb=1.0, beta=1.0):
	F=1.0-scipy.array(map(math.exp, -(x/lamb)**beta))
	#zs=1.0-scipy.array(map(math.exp, -(scipy.array(ints)/tbar)**beta))
	return F

def randomDists(N=10000):
	r1=random.Random()
	Rs=[]
	for i in xrange(N):
		Rs+=[r1.random()]
	Rs=scipy.array(Rs)
	
	Rsums = (Rs[0:-1] + Rs[1:])
	Rdifs = (Rs[0:-1] - Rs[1:])
	
	for x in range(4):
		plt.figure(x)
		plt.clf()
	#return None
	#
	'''
	plt.figure(0)
	plt.plot(Rsums, '.')
	plt.figure(1)
	plt.plot(Rdifs, '.')
	'''
	#
	plt.figure(0)
	plt.hist(Rsums, bins=100, normed=True)
	plt.figure(1)
	plt.hist(Rdifs, bins=100, normed=True)
	
def rbmarkov(catname='cats/socaliperc.cat', lon=[-122.5, -114.0], lat=[30.5, 37.25], minMag=2.5, dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None], Nmax=999999, beta1=.5, beta2=1.0, refreshcat=True, dataset=None, winlen=320):
	c1=getCat(catname=catname, lon=lon, lat=lat, minMag=minMag, dates0=dates0, Nmax=Nmax, beta1=beta1, beta2=beta2, refreshcat=refreshcat, dataset=dataset)
	rbints=c1.getNRBratios(winlen=winlen)
	avelen=int(winlen/10)
	#
	logZ=math.log10(winlen)
	#
	r_vals=map(operator.itemgetter(4), rbints)
	log_rvals=scipy.array(map(math.log10, r_vals))
	log_rvals/=logZ
	delta_rs=log_rvals[1:]-log_rvals[0:-1]
	fdts=map(mpd.date2num, map(operator.itemgetter(1), rbints))
	#
	r_vals_averaged=scipy.array(averageSequence(log_rvals, avelen=avelen))
	delta_rsAveraged = r_vals_averaged[1:]-r_vals_averaged[:-1]
	#
	f0=plt.figure(0)
	plt.clf()
	axdn=f0.add_axes([.05, .05, .85, .45])
	axup=f0.add_axes([.05, .5, .85, .45], sharex=axdn)
	axup.plot(fdts, log_rvals, '.', zorder=1)
	axdn.plot(fdts[1:], delta_rs, '.', zorder=1)
	axup.plot(fdts[avelen:], r_vals_averaged, '.', zorder=2, alpha=.5)
	axdn.plot(fdts[1+avelen:], delta_rsAveraged, '.', zorder=2, alpha=.5)
	#
	f0=plt.figure(1)
	plt.clf()
	axdn=f0.add_axes([.05, .05, .85, .45])
	axup=f0.add_axes([.05, .5, .85, .45], sharex=axdn)
	axup.plot(log_rvals, '.', zorder=1)
	axdn.plot(delta_rs, '.', zorder=1)
	axup.plot(r_vals_averaged, '.', zorder=2, alpha=.5)
	axdn.plot(delta_rsAveraged, '.', zorder=2, alpha=.5)
	
	#
	f1=plt.figure(2)
	plt.clf()
	axleft=f1.add_axes([.05, .1, .45, .85])
	axright=f1.add_axes([.5, .1, .45, .85])
	axright2=axright.twinx()
	h1=axleft.hist(log_rvals, bins=100, normed=True)
	axleft.set_xlabel('log_rvals')
	h2=axright.hist(delta_rs, bins=100, normed=True)
	axright.set_xlabel('$\\Delta [\\log(r)]$')
	ldrcdf=delta_rs.copy()
	ldrcdf.sort()
	Ns=range(1, len(ldrcdf)+1)
	axright2.plot(ldrcdf, Ns, 'r-')
	#
	axleft.plot([0., 0.], [0., 1.2*max(h1[0])], 'm-',lw=2, alpha=.7)
	axright.plot([0., 0.], [0., 1.2*max(h2[0])], 'm-',lw=2, alpha=.7)
	
	f1=plt.figure(3)
	plt.clf()
	axleft=f1.add_axes([.05, .1, .45, .85])
	axright=f1.add_axes([.5, .1, .45, .85])
	axright2=axright.twinx()
	#
	h1=axleft.hist(r_vals_averaged, bins=100, normed=True)
	axleft.set_xlabel('log_rvals averaged')
	h2=axright.hist(delta_rsAveraged, bins=100, normed=True)
	axright.set_xlabel('$\\Delta [\\log(r)]$')
	ldrcdf=delta_rsAveraged.copy()
	ldrcdf.sort()
	Ns=range(1, len(ldrcdf)+1)
	axright2.plot(ldrcdf, Ns, 'r-')
	#
	axleft.plot([0., 0.], [0., 1.2*max(h1[0])], 'm-',lw=2, alpha=.7)
	axright.plot([0., 0.], [0., 1.2*max(h2[0])], 'm-',lw=2, alpha=.7)
	#
	# treat r>0 and r<0 separately?
	rgt=[]
	rlt=[]
	for rval in log_rvals:
		if rval>0.: rgt+=[rval]
		if rval<0.: rlt+=[rval]
	rgt=scipy.array(rgt)
	rlt=scipy.array(rlt)
	deltargt=rgt[1:]-rgt[:-1]
	deltarlt=rlt[1:]-rlt[:-1]
	#
	f1=plt.figure(4)
	plt.clf()
	axll=f1.add_axes([.05, .1, .45, .45])
	axlr=f1.add_axes([.5, .1, .45, .45])
	axul=f1.add_axes([.05, .5, .45, .45])
	axur=f1.add_axes([.5, .5, .45, .45])
	#axright2=axright.twinx()
	h1=axll.hist(rgt, bins=100, normed=True)
	h2=axlr.hist(deltargt, bins=100, normed=True)
	h3=axul.hist(rlt, bins=100, normed=True)
	h4=axur.hist(deltarlt, bins=100, normed=True)
	#
	axll.set_xlabel('log_rvals gt')
	axlr.set_xlabel('delta rs, gt')
	axul.set_xlabel('log_rvals lt')
	axur.set_xlabel('delta rs, lt')
	
	ldrcdf=delta_rs.copy()
	ldrcdf.sort()
	Ns=range(1, len(ldrcdf)+1)
	axright2.plot(ldrcdf, Ns, 'r-')
	#
	axleft.plot([0., 0.], [0., 1.2*max(h1[0])], 'm-',lw=2, alpha=.7)
	axright.plot([0., 0.], [0., 1.2*max(h2[0])], 'm-',lw=2, alpha=.7)


def averageSequence(X, avelen=1):
	aveX=[]
	i=avelen
	while i<len(X):
		aveX+=[scipy.mean(X[i-avelen:i])]
		i+=1
	return aveX

def getCat(catname='cats/socaliperc.cat', lon=[-122.5, -114.0], lat=[30.5, 37.25], minMag=2.5, dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None], Nmax=999999, beta1=.5, beta2=1.0, refreshcat=True, dataset=None):
	# catfromANSS(lon=[135., 150.], lat=[30., 41.5], minMag=4.0, dates0=[dtm.datetime(2005,01,01, tzinfo=tzutc), None], Nmax=999999, fout='cats/mycat.cat')
	#tohoku:
	# lons=[135., 146.], lats=[30., 41.5]
	# approx for norcal:
	# ipp.intervalDistAnalysis(catname='cats/norcaliperc.cat', lon=[-128., -117.], lat=[37., 43.], minMag=3.0, refreshcat=True)
	if dataset=='tohoku':
		catname='cats/tohokutmp.cat'
		lon=[135., 150.]
		lat=[30., 41.5]
		minMag=4.5
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), dtm.datetime(2011,03,01, tzinfo=pytz.timezone('UTC'))]
		dates0=[dtm.datetime(2011,3,11, 12, tzinfo=pytz.timezone('UTC')), None]
	if dataset=='socal':
		catname='cats/socaliperc.cat'
		lon=[-122.5, -114.0]
		lat=[30.5, 37.25]
		minMag=2.5
		#dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
		dates0=[dtm.datetime(2009,01,01, tzinfo=pytz.timezone('UTC')), dtm.datetime(2010,01,01, tzinfo=pytz.timezone('UTC'))]
	if dataset=='norcal':
		catname='cats/norcaliperc.cat'
		lon=[-128., -117.]
		lat=[37., 43.]
		minMag=3.0
		dates0=[dtm.datetime(2000,01,01, tzinfo=pytz.timezone('UTC')), None]
	#
	if refreshcat:
		clist = atp.catfromANSS(lon=lon, lat=lat, minMag=minMag, dates0=dates0, Nmax=Nmax, fout=catname)
		clist.sort(key=lambda x:x[0])
		c1=eqp.eqcatalog(clist)
	else:
		c1=eqp.eqcatalog()
		c1.loadCatFromFile(catname)
	#
	return c1	
	
