"""Compute maximum entropy coefficients for data. Based loosely on the the
concepts in "Numerical Recipes in C", 2nd ed., by Press, Flannery, Teukolsky and Vetterling (q.v.)
but I don't think it is any copyright infringement"""
_rcsid="$Id: memcof.py,v 1.3 2003-05-30 13:31:55 mendenhall Exp $"

import Numeric
import math

def memcof(data, poles):
	n=len(data)
	xms=Numeric.dot(data,data)/float(n)
	wk1=data[:-1]
	wk2=data[1:]
	d=Numeric.zeros(poles,Numeric.Float)
	
	for k0 in range(poles):
		num=Numeric.dot(wk1,wk2)
		denom=Numeric.dot(wk1, wk1)+Numeric.dot(wk2, wk2)
		d[k0]=2.0*num/denom
		xms*=(1.0-d[k0]**2)
		if k0!=0:
			d[:k0]=wkm-d[k0]*wkm[-1::-1]			
		if k0!=poles-1:
			wkm=d[:k0+1]
			wk1, wk2 = wk1[:-1]-wkm[k0]*wk2[:-1], wk2[1:]-wkm[k0]*wk1[1:]

	return xms, d


def evlmem(fdt, d, xms):
	n=len(d)
	theta=2*math.pi*fdt*(Numeric.array(range(1,n+1), Numeric.Float))
	zr=1.0-Numeric.dot(Numeric.cos(theta), d)
	zi=Numeric.dot(Numeric.sin(theta), d)
	return xms/(zr*zr+zi*zi)
	

if __name__=="__main__":
	
	datalen=10000
	poles=20
	zfreq=0.21345
	pspoints=200
	damping=-20.0/datalen
	noise=1.0
	
	xvals=Numeric.array(range(datalen),Numeric.Float)
	data=Numeric.sin(xvals*(2.0*math.pi*zfreq)) * Numeric.exp(xvals*damping)
	import random
	r=random.Random(1234)
	data+= Numeric.array([r.random()-0.5 for i in range(datalen)])*noise
	
	d2=Numeric.dot(data,data)/datalen #mean-square power
	
	xms, d = memcof(data, poles)
	
	freqlist = [0.5*i/pspoints for i in range(pspoints)]
	pspect = [evlmem(f, d, xms) for f in freqlist]
	
	pssum=2.0*Numeric.sum(pspect)*(0.5/pspoints)
	
	print "input power = ", d2, "output power = ", pssum, "ratio =", pssum/d2
	import graphite
	g=graphite.Graph()
	d=graphite.Dataset()
	d.x=freqlist
	d.y=pspect
	g.datasets=d
	graphite.genOutput(g,'QD',canvasname="Power Spectrum")
	
