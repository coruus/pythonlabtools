"""cubic spline handling, in a manner compatible with the API in Numeric Recipes"""
_rcsid="$Id: spline.py,v 1.8 2003-06-26 18:19:06 mendenhall Exp $"

__all__=["spline","splint","cubeinterpolate","RangeError"]

class RangeError(IndexError):
	"X out of input range in splint()"

from Numeric import zeros, Float, searchsorted, array, asarray, take, clip

def spline(x, y, yp1=None, ypn=None):
	"""spline(x_vals,y_vals, yp1=None, ypn=None) 
	returns the y2 table for the spline as needed by splint()"""

	n=len(x)
	u=zeros(n,Float)
	y2=zeros(n,Float)
	
	x=asarray(x, Float)
	y=asarray(y, Float)
	
	dx=x[1:]-x[:-1]
	dxi=1.0/dx
	dx2i=1.0/(x[2:]-x[:-2])
	dy=(y[1:]-y[:-1])
	siga=dx[:-1]*dx2i
	dydx=dy*dxi
	
	# u[i]=(y[i+1]-y[i])/float(x[i+1]-x[i]) - (y[i]-y[i-1])/float(x[i]-x[i-1])
	u[1:-1]=dydx[1:]-dydx[:-1] #this is an incomplete rendering of u... the rest requires recursion in the loop
	
	if yp1 is None:
		y2[0]=u[0]=0.0
	else:
		y2[0]= -0.5
		u[0]=(3.0*dxi[0])*(dy[0]*dxi[0] -yp1)

	for i in range(1,n-1):
		sig=siga[i-1]
		p=sig*y2[i-1]+2.0
		y2[i]=(sig-1.0)/p
		u[i]=(6.0*u[i]*dx2i[i-1] - sig*u[i-1])/p

	if ypn is None:
		qn=un=0.0
	else:
		qn= 0.5
		un=(3.0*dxi[-1])*(ypn- dy[-1]*dxi[-1] )
		
	y2[-1]=(un-qn*u[-2])/(qn*y2[-2]+1.0)
	for k in range(n-2,-1,-1):
		y2[k]=y2[k]*y2[k+1]+u[k]

	return y2

import types

def splint(xa, ya, y2a, x):
	"""splint(x_vals, y_vals, y2_vals, x) returns the interpolated from from the spline
	x can either be a scalar or a listable item, in which case a Numeric Float array will be
	returned and the multiple interpolations will be done somewhat more efficiently if the
	array is in reasonable order"""
	if type(x) is types.IntType or type(x) is types.FloatType: 
		if x<xa[0] or x>xa[-1]:
			raise RangeError, "%f not in range (%f, %f) in splint()" % (x, xa[0], xa[-1])
			 
		khi=max(searchsorted(xa,x),1)
		klo=khi-1
		h=float(xa[khi]-xa[klo])
		a=(xa[khi]-x)/h; b=1.0-a
		return a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0
	else:
		#if we got here, we are processing a list, and should do so more efficiently
		if min(x)<xa[0] or max(x)>xa[-1]:
			raise RangeError, "(%f, %f) not in range (%f, %f) in splint()" % (min(x), max(x), xa[0], xa[-1])
	
		npoints=len(x)

		khi=clip(searchsorted(xa,x),1,npoints) 
		klo=khi-1
		xhi=take(xa, khi)
		xlo=take(xa, klo)
		yhi=take(ya, khi)
		ylo=take(ya, klo)
		y2hi=take(y2a, khi)
		y2lo=take(y2a, klo)
		
		h=(xhi-xlo).astype(Float)
		a=(xhi-x)/h
		b=1.0-a
		
		return a*ylo+b*yhi+((a*a*a-a)*y2lo +(b*b*b-b)*y2hi)*(h*h)/6.0
		
def cubeinterpolate(xlist, ylist, x3):
	"find point at x3 given 4 points in given lists using exact cubic interpolation, not splining"
	x1,x2,x4,x5=xlist
	x2,x3,x4,x5=float(x2-x1),float(x3-x1),float(x4-x1),float(x5-x1)
	y1,y2,y4,y5=ylist
	y2,y4, y5=float(y2-y1),float(y4-y1),float(y5-y1)
	
	y3=(
			(x3*(x2**2*x5**2*(-x2 + x5)*y4 + x4**3*(x5**2*y2 - x2**2*y5) + x4**2*(-(x5**3*y2) + x2**3*y5) + 
		           x3**2*(x2*x5*(-x2 + x5)*y4 + x4**2*(x5*y2 - x2*y5) + x4*(-(x5**2*y2) + x2**2*y5)) + 
		           x3*(x2*x5*(x2**2 - x5**2)*y4 + x4**3*(-(x5*y2) + x2*y5) + x4*(x5**3*y2 - x2**3*y5))))/
		      	 (x2*(x2 - x4)*x4*(x2 - x5)*(x4 - x5)*x5)
	)+y1
	return y3


if __name__=="__main__":
	import traceback
	testlist=((0,1), (1,1),(2,3),(3,4),(4,2),(5,6),(7,9),(10,6),(15,2), (16,-1))
	#testlist=((0,0), (1,1),(2,4),(3,9),(4,16),(5,25),(7,49),(10,100),(15,225), (16,256))
	xlist=[i[0] for i in testlist]
	ylist=[i[1] for i in testlist]
	print "\n\nStarting splint tests...\n", testlist
	y2=spline(xlist,ylist, yp1=-5, ypn=10)
	r=(0,1,2,3.5, 3.7, 4,6,7,2,8,9,10,11, 5, 12,13,14, 15)
	v=splint(xlist, ylist, y2, r)	
	print y2
	for i in range(len(r)):
		print "%.1f %.3f %.3f" % (r[i], v[i], splint(xlist, ylist, y2, r[i]))
	
	print "The next operations should print exceptions"
	try:
		splint(xlist, ylist, y2, 100.0)
	except:
		traceback.print_exc()
	try:
		splint(xlist, ylist, y2, (1,2,2.5, 3,-5, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18))
	except:
		traceback.print_exc()
	
	try:
		import graphite
		import Numeric
		g=graphite.Graph()
		ds1=graphite.Dataset()
		ds1.x=[i[0] for i in testlist]
		ds1.y=[i[1] for i in testlist]
		g.datasets.append(ds1)
		f1 = graphite.PointPlot()
		f1.lineStyle = None
		f1.symbol = graphite.CircleSymbol
		f1.symbolStyle=graphite.SymbolStyle(size=5, fillColor=graphite.red, edgeColor=graphite.red)
		g.formats=[]
		g.formats.append(f1)
		finex=Numeric.array(range(0,161),Float)*0.1
		finey=splint(xlist, ylist, y2, finex)
		ds2=graphite.Dataset()
		ds2.x=finex
		ds2.y=finey
		g.datasets.append(ds2)
		f2 = graphite.PointPlot()
		f2.lineStyle = graphite.LineStyle(width=1, color=graphite.green, kind=graphite.SOLID)
		f2.symbol = None
		g.formats.append(f2)
		g.bottom=400
		g.right=700
		try:
			graphite.genOutput(g,'QD', size=(800,500))
		except:
			graphite.genOutput(g,'PDF', size=(800,500))
	except:
		import traceback
		traceback.print_exc()
		print "Graphite not available... plotted results not shown"
		
