"""cubic spline handling, in a manner compatible with the API in Numeric Recipes"""
#$Id: spline.py,v 1.2 2003-04-11 18:06:16 mendenhall Exp $

import exceptions

class RangeError(exceptions.IndexError):
	"X out of input range in splint()"

from Numeric import zeros, Float, searchsorted, array

def spline(x, y, yp1=None, ypn=None):
	"""spline(x_vals,y_vals, yp1=None, ypn=None) 
	returns the y2 table for the spline as needed by splint()"""

	n=len(x)
	u=zeros(n,Float)
	y2=zeros(n,Float)
	
	
	if yp1 is None:
		y2[0]=u[0]=0.0
	else:
		y2[0]= -0.5
		u[0]=(3.0/(x[1]-x[0]))*( ((y[1]-y[0])/float(x[1]-x[0])) -yp1)

	for i in range(1,n-1):
		sig=(x[i]-x[i-1])/float(x[i+1]-x[i-1])
		p=sig*y2[i-1]+2.0
		y2[i]=(sig-1.0)/p
		u[i]=(y[i+1]-y[i])/float(x[i+1]-x[i]) - (y[i]-y[i-1])/float(x[i]-x[i-1])
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1])/p

	if ypn is None:
		qn=un=0.0
	else:
		qn= 0.5
		un=(3.0/(x[-1]-x[-2]))*(ypn- ((y[-1]-y[-2])/float(x[-1]-x[-2])) )
		
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
		results=zeros(npoints, Float)
		n=len(xa)
		klo=0
		khi=1
		for i in range(len(x)):
			xv=x[i]
			if not (xa[klo] <= xv < xa[klo+1]): #number is not currently in range
				bound=min(klo+5, n)		
		 		khi=klo+searchsorted(xa[klo:bound], xv) #first try local search for efficiency
				if khi==bound or khi==klo:
					khi=max(searchsorted(xa, xv),1) #local didn't work, try whole array so big jumps / out of order works
				klo=khi-1
			h=float(xa[khi]-xa[klo])
			a=(xa[khi]-xv)/h; b=1.0-a
			results[i]=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0
		return results

if __name__=="__main__":
	import traceback
	testlist=((0,1), (1,1),(2,3),(3,4),(4,2),(5,6),(7,9),(10,6),(15,2), (16,-1))
	#testlist=((0,0), (1,1),(2,4),(3,9),(4,16),(5,25),(7,49),(10,100),(15,225), (16,256))
	xlist=[i[0] for i in testlist]
	ylist=[i[1] for i in testlist]
	print "\n\nStarting splint tests...\n", testlist
	y2=spline(xlist,ylist, yp1=10, ypn=0)
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
			graphite.genOutput(q,'PDF', size=(800,500))
	except:
		print "Graphite not available... plotted results not shown"
		