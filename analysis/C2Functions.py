"""A group of classes which make it easy to manipulate smooth functions, including cubic splines. 

C2Functions know how to keep track of the first and second derivatives of functions, and to use this iformation in, for example, find_root()
to allow much more efficient solutions to problems for which the general solution may be expensive.

The two primary classes are 
	C2Function, which represents an unevaluted function and its derivatives, and 
	InterpolatingFunction,  which represent a cubic spline of provided data.

C2Functions can be combined with unary operators (nested functions) or binary operators (+-*/ etc.)

Developed by Marcus H. Mendenhall, Vanderbilt University Keck Free Electron Laser Center, Nashville, TN USA
email: marcus.h.mendenhall@vanderbilt.edu
Work supported by the US DoD  MFEL program under grant FA9550-04-1-0045
"""
_rcsid="$Id: C2Functions.py,v 1.1 2005-07-18 17:00:33 mendenhall Exp $"

import math
import operator
import types

import Numeric as _numeric #makes it easy to change to NumArray later
_myfuncs=_numeric #can change to math for (possibly greater) speed (?) but no vectorized evaluation

class C2Exception(Exception):
	pass
	
class C2Function:
	"if f is a C2Function, f(x) returns the value at x, and f.value_with_derivatives returns y(x), y'(x), y''(x)"
	ClassName='C2Function'
	name='Unnamed'
	def __init__(self, *args) : 
		if not args:
			self.xMin, self.xMax=-1e308, 1e308
		elif len(args)==1: #constructing from a single C2Function
			self.xMin, self.xMax=args[0].xMin, args[0].xMax
		else: #it should be two C2Functions, and we take the inner bounds of the two
			l,r=args
			self.xMin=max(l.xMin, r.xMin)
			self.xMax=min(l.xMax, r.xMax)
		
		
	def __str__(self): return '<%s %s, Domain=[%g, %g]>' % ((self.ClassName, self.name,)+ self.GetDomain())
	
	def SetName(self, name): 
		"set the short name of the function for __str__"
		self.name=name; return self
	
	def __call__(self, arg): 
		"f(x) evaluates f(arg) without derivatives if arg is numeric, otherwise returns the composed function f(arg())"
		if isinstance(arg, C2Function): return C2ComposedFunction(self, arg)
		else: return self.value_with_derivatives(arg)[0] 

	def apply(self, interpolator):
		"creates a new InterpolatingFunction which has this function applied to the original interpolater"
		return interpolator.UnaryOperator(self)

	def find_root(self, lower_bracket, upper_bracket, start, value, trace=False): # solve f(x)=value
		"solve f(x)==value very efficiently, with explicit knowledge of derivatives of the function"
		#find f(x)=value within the brackets, using the guarantees of smoothness associated with a C2Function
		# can use local f(x)=a*x**2 + b*x + c and solve quadratic to find root, then iterate

		ftol=(2e-16*abs(value)+1e-308);
		xtol=(2e-16*(abs(upper_bracket)+abs(lower_bracket))+1e-308);

		root=start	#start looking in the middle
		cupper, yp, ypp=self.value_with_derivatives(upper_bracket)
		cupper -= value
		if abs(cupper) < ftol: return upper_bracket
		clower, yp, ypp=self.value_with_derivatives(lower_bracket)
		clower -= value
		if abs(clower) < ftol: return lower_bracket;

		if cupper*clower >0:
			# argh, no sign change in here!
			raise C2Exception("unbracketed root in find_root at xlower=%g, xupper=%g, bailing" %(lower_bracket, upper_bracket))

		delta=upper_bracket-lower_bracket	#first error step
		c, b, ypp=self.value_with_derivatives(root) # compute initial values
		c -= value
		
		while abs(delta) > xtol: # can allow quite small steps, since we have exact derivatives!
			a=ypp/2	#second derivative is 2*a
			disc=b*b-4*a*c
			if disc > 0:
				if b>=0:
					q=-0.5*(b+math.sqrt(disc))
				else:
					q=-0.5*(b-math.sqrt(disc))

				if q*q > abs(a*c): delta=c/q	#since x1=q/a, x2=c/q, x1/x2=q^2/ac, this picks smaller step
				else: delta=q/a
				root+=delta;

			if disc <= 0 or root<lower_bracket or root>upper_bracket:	#if we jump out of the bracket, bisect
				root=0.5*(lower_bracket+upper_bracket)	
				delta=upper_bracket-lower_bracket

			c, b, ypp=self.value_with_derivatives(root)	#get local curvature and value
			c -= value

			if trace: 
				import sys
				print >> sys.stderr, "find_root x, dx, c, b, a", (5*"%10.4g ") % (root, delta, c, b, ypp/2)
				
			if abs(c) < ftol or abs(c) < xtol*abs(b):	return root		#got it close enough

			# now, close in bracket on whichever side this still brackets
			if c*clower < 0.0:
				cupper=c
				upper_bracket=root
			else:
				clower=c
				lower_bracket=root

		return root

	def SetDomain(self, xmin, xmax): 
		"Set the domain of the function. This is mostly an advisory range, except for InterpolatingFunctions, where it matters"
		self.xMin=xmin; self.xMax=xmax; return self
	def GetDomain(self):
		"returns xMin, xMax"
		return self.xMin, self.xMax

	def compose(self, inner_function, x):
		"y=f.compose(g, x) returns f(g(x)), f'(g(x)), f''(g(x)) where the derivatives are with respect to x"
		y0, yp0, ypp0=inner_function.value_with_derivatives(x)
		y1, yp1, ypp1=self.value_with_derivatives(y0)
		return y1, yp1*yp0, ypp0*yp1+yp0*yp0*ypp1

	def convert_arg(self, arg): 
		"convert constants to C2Constants"
		import operator
		if type(arg) is types.FloatType or type(arg) is types.IntType:
			return C2Constant(arg)
		else: 
			return arg

	def __add__(self, right):
		"a+b returns a new C2Function which represents the sum"
		return C2Sum(self, right)
	def __sub__(self, right):
		"a-b returns a new C2Function which represents the difference"
		return C2Diff(self, right)
	def __mul__(self, right):
		"a*b returns a new C2Function which represents the product"
		return C2Product(self, right)
	def __div__(self, right):
		"a/b returns a new C2Function which represents the ratio"
		return C2Ratio(self, right)

class C2Constant(C2Function):
	"a constant and its derivatives"
	ClassName='Constant'
	def __init__(self, val):
		C2Function.__init__(self)
		self.val=val
		self.name='%g' % val
	def value_with_derivatives(self, x): 
		return self.val, 0., 0.

class _fC2sin(C2Function):
	"sin(x)"
	name='sin'
	def value_with_derivatives(self, x):
		q=_myfuncs.sin(x)
		return q, _myfuncs.cos(x), -q
C2sin=_fC2sin() #make singleton

class _fC2cos(C2Function):
	"cos(x)"
	name='cos'
	def value_with_derivatives(self, x):
		q=_myfuncs.cos(x)
		return q, -_myfuncs.sin(x), -q
C2cos=_fC2cos() #make singleton

class _fC2log(C2Function):
	"log(x)"
	name='log'
	def value_with_derivatives(self, x):
		return _myfuncs.log(x), 1/x, -1/(x*x)
C2log=_fC2log() #make singleton

class _fC2exp(C2Function):
	"exp(x)"
	name='exp'
	def value_with_derivatives(self, x):
		q=_myfuncs.exp(x)
		return q, q, q
C2exp=_fC2exp() #make singleton

class _fC2sqrt(C2Function):
	"sqrt(x)"
	name='sqrt'
	def value_with_derivatives(self, x):
		q=_myfuncs.sqrt(x)
		return q, 0.5/q, -0.25/(q*x)
C2sqrt=_fC2sqrt() #make singleton

class _fC2recip(C2Function):
	"1/x"
	name='1/x'
	def value_with_derivatives(self, x):
		q=1.0/x
		return q, -q*q, 2*q*q*q
C2recip=_fC2recip() #make singleton

class _fC2identity(C2Function):
	name='Identity'
	def value_with_derivatives(self, x):
		return x, 1.0, 0.0
C2identity=_fC2identity() #make singleton

class C2Linear(C2Function):
	"slope*x + y0"
	def __init__(self, slope=1.0, y0=0.0):
		C2Function.__init__(self)
		self.slope=slope
		self.y0=y0
		self.name="(%g * x + %g)" % (slope, y0)
		
	def value_with_derivatives(self, x):
		return x*self.slope+self.y0, self.slope, 0.0

class C2Quadratic(C2Function):
	"a*(x-x0)**2 + b*(x-x0) + c"
	def __init__(self, x0=0.0, a=1.0, b=0.0, c=0.0):
		C2Function.__init__(self)
		self.x0, self.a, self.b, self.c = x0, a, b, c
		self.name="(%g*(x-x0)^2 + %g*(x-x0) + %g, x0=%g)" % (a, b, c, x0)
		
	def value_with_derivatives(self, x):
		dx=x-self.x0
		return self.a*dx*dx+self.b*dx+self.c, 2*self.a*dx+self.b, 2*self.a

class C2PowerLaw(C2Function):
	"a*x**b"
	def __init__(self, a=1.0, b=2.0):
		C2Function.__init__(self)
		self.a, self.b = a, b
		self.b2=b-2
		self.bb1=b*(b-1)
		self.name='%g*x^%g' % (a,b)
		
	def value_with_derivatives(self, x):
		xp=self.a*x**self.b2
		return xp*x*x, self.b*xp*x, self.bb1*xp

class C2ComposedFunction(C2Function):
	"create a composed function outer(inner(...)).  The functions can either be unbound class names or class instances"
	def __init__(self, outer, inner):
		if type(inner) is types.ClassType: inner=inner() #instatiate unbound class
		if type(outer) is types.ClassType: outer=outer() #instatiate unbound class
		C2Function.__init__(self, inner) #domain is that of inner function
		self.outer=outer
		self.inner=inner
		self.name=outer.name+'('+inner.name+')'
		
	def value_with_derivatives(self, x): return self.outer.compose(self.inner, x)

class C2BinaryFunction(C2Function):
	"create a binary function using a helper function which computes the derivatives"
	def __init__(self, left,  right):
		C2Function.__init__(self, left, right)
		self.left=self.convert_arg(left)
		self.right=self.convert_arg(right)
		
		if isinstance(left, C2BinaryFunction): p1, p2 = '(', ')'
		else: p1, p2='', ''
		self.name=p1+left.name+p2+self.name+right.name #put on parentheses to kepp hierachy obvious
		
	def value_with_derivatives(self, x):
		return self.combine(x);

class C2Sum(C2BinaryFunction):
	"C2Sum(a,b) returns a new C2Function which evaluates as a+b"
	name='+'
	def combine(self, x):
		y0, yp0, ypp0=self.left.value_with_derivatives(x)
		y1, yp1, ypp1=self.right.value_with_derivatives(x)
		return y0+y1, yp0+yp1, ypp0+ypp1

class C2Diff(C2BinaryFunction):
	"C2Diff(a,b) returns a new C2Function which evaluates as a-b"
	name='-'
	def combine(self, x):
		y0, yp0, ypp0=self.left.value_with_derivatives(x)
		y1, yp1, ypp1=self.right.value_with_derivatives(x)
		return y0-y1, yp0-yp1, ypp0-ypp1

class C2Product(C2BinaryFunction):
	"C2Product(a,b) returns a new C2Function which evaluates as a*b"
	name='*'
	def combine(self, x):
		y0, yp0, ypp0=self.left.value_with_derivatives(x)
		y1, yp1, ypp1=self.right.value_with_derivatives(x)
		return y0*y1, y1*yp0+y0*yp1, ypp0*y1+2.0*yp0*yp1+ypp1*y0

class C2Ratio(C2BinaryFunction):
	"C2Ratio(a,b) returns a new C2Function which evaluates as a/b"
	name='/'
	def combine(self, x):
		y0, yp0, ypp0=self.left.value_with_derivatives(x)
		y1, yp1, ypp1=self.right.value_with_derivatives(x)
		return y0/y1, (yp0*y1-y0*yp1)/(y1*y1), (y1*y1*ypp0+y0*(2*yp1*yp1-y1*ypp1)-2*y1*yp0*yp1)/(y1*y1*y1)

class InterpolatingFunction(C2Function):
	"""An InterpolatingFunction stores a cubic spline representation of a set of x,y pairs.
		It can also transform the variable on input and output, so that the underlying spline may live in log-log space, 
		but such transforms are transparent to the setup and use of the function.  This makes it possible to
		store splines of, e.g., data which are very close to a power law, as a LogLogInterpolatingFunction, and
		to then have very accurate interpolation and extrapolation, since the curvature of such a function is small in log-log space.
		
		InterpolatingFunction(x, y, lowerSlope, upperSlope, XConversions, YConversions) sets up a spline.  
		If lowerSlope or upperSlope is None, the corresponding boundary is set to 'natural', with zero second derivative.
		XConversions is a list of g, g', g'' to evaluate for transforming the X axis.  
		YConversions is a list of f, f', f'', f(-1) to evaluate for transforming the Y axis.
			Note that the y transform f and f(-1) MUST be exact inverses, or the system will melt.
				
	""" 
	YConversions=None
	XConversions=None
	name='data'
	ClassName='InterpolatingFunction'
	def __init__(self, x, y, lowerSlope=None, upperSlope=None, XConversions=None, YConversions=None):
		C2Function.__init__(self) #just on general principle, right now this does nothing
		self.SetDomain(min(x), max(x))
		self.Xraw=_numeric.array(x) #keep a private copy

		if YConversions is not None: self.YConversions=YConversions #inherit from class if not passed		
		if self.YConversions is None:
			self.fYin, self.fYinP, self.fYinPP, self.fYout = self.YConversions = lambda x: x, lambda x: 1, lambda x: 0, lambda x : x
			self.yNonLin=False
			F=self.F=_numeric.array(y)
		else:
			self.fYin, self.fYinP, self.fYinPP, self.fYout = self.YConversions
			self.yNonLin=True
			self.F=_numeric.array([self.fYin(q) for q in y])
			if lowerSlope is not None: lowerSlope *= self.fYinP(y[0])
			if upperSlope is not None: upperSlope *= self.fYinP(y[-1])

		if XConversions is not None: self.XConversions=XConversions #inherit from class if not passed		
		if self.XConversions is None:
			self.fXin, self.fXinP, self.XinPP = self.XConversions =  lambda x: x, lambda x: 1, lambda x: 0
			self.xNonLin=False
			self.X=_numeric.array(x)
		else:
			self.fXin, self.fXinP, self.fXinPP = self.XConversions
			self.xNonLin=True
			self.X=_numeric.array([self.fXin(q) for q in x])
			if lowerSlope is not None: lowerSlope /= self.fXinP(x[0])
			if upperSlope is not None: upperSlope /= self.fXinP(x[-1])

			if self.X[0] > self.X[-1]: #make sure our transformed X array is increasing
				self.Xraw=self.Xraw[::-1]
				self.X=self.X[::-1]
				self.Y=self.Y[::-1]
				lowerSlope, upperSlope=upperSlope, lowerSlope

		from analysis import spline
		self.y2=spline.spline(self.X, self.F, yp1=lowerSlope, ypn=upperSlope)

	def value_with_derivatives(self, x): 
		from analysis import spline
		if self.xNonLin or self.yNonLin: #skip this if this is a completely untransformed spline
			y0, yp0, ypp0=spline.splint(self.X, self.F, self.y2, self.fXin(x), derivs=True)
			y=self.fYout(y0)
			fpi=1.0/self.fYinP(y)
			gp=self.fXinP(x)

			# from Mathematica Dt[InverseFunction[f][g[y[x]]], x]
			yprime=gp*yp0*fpi # the real derivative of the inverse transformed output 
			fpp=self.fYinPP(y)
			gpp=self.fXinPP(x);
			#also from Mathematica Dt[InverseFunction[f][g[y[x]]], {x,2}]
			yprimeprime=(gp*gp*ypp0 + yp0*gpp - gp*gp*yp0*yp0*fpp*fpi*fpi)*fpi; 
			return y, yprime, yprimeprime
		else:
			return spline.splint(self.X, self.F, self.y2, x, derivs=True)

	def UnaryOperator(self, C2source):
		"return new InterpolatingFunction C2source(self)"
		yv=[C2source.compose(self, x) for x in self.Xraw] #get array of ( ( y, y', y''), ...)
		y=[yy[0] for yy in yv] #get y values only

		f=InterpolatingFunction(self.Xraw, y, lowerSlope=yv[0][1], upperSlope=yv[-1][1],
			XConversions=self.XConversions, YConversions=self.YConversions)
		f.name=C2source.name+'('+self.name+')'
		return f
		
	def BinaryOperator(self, rhs, c2binary):
		"return new InterpolatingFunction self +-*/ rhs (or any other binary operator)"
		bf=c2binary(self, rhs)

		yv=[bf.value_with_derivatives(x) for x in self.Xraw] #get array of ( ( y, y', y''), ...)
		y=[yy[0] for yy in yv] #get y values only

		f=InterpolatingFunction(self.Xraw, y, lowerSlope=yv[0][1], upperSlope=yv[-1][1],
			XConversions=self.XConversions, YConversions=self.YConversions)
		f.name=bf.name
		return f
		
	#factory functions to create evaluated binary functions of InterpolatingFunctions
	def __add__(self, right):
		return self.BinaryOperator(right, C2Sum)
	def __sub__(self, right):
		return self.BinaryOperator(right, C2Diff)
	def __mul__(self, right):
		return self.BinaryOperator(right, C2Product)
	def __div__(self, right):
		return self.BinaryOperator(right, C2Ratio)

class LogLinInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores log(x) vs. y"
	ClassName='LogLinInterpolatingFunction'
	XConversions=math.log, lambda x: 1.0/x, lambda x: -1.0/(x*x)

class LinLogInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores x vs. log(y), useful for functions with exponential-like behavior"
	ClassName='LinLogInterpolatingFunction'
	YConversions=math.log, lambda x: 1.0/x, lambda x: -1.0/(x*x), math.exp

class LogLogInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores log(x) vs. log(y), useful for functions with power-law-like behavior"
	ClassName='LogLogInterpolatingFunction'
	XConversions=math.log, lambda x: 1.0/x, lambda x: -1.0/(x*x)
	YConversions=math.log, lambda x: 1.0/x, lambda x: -1.0/(x*x), math.exp

def LinearInterpolatingGrid(xmin, dx, count):
	"""create a linear-linear interpolating grid with both x & y set to (xmin, xmin+dx, ... xmin + dx*(count -1) )
		very useful for transformaiton with other functions e.g. 
		f=C2sin(LinearInterpolatingGrid(-0.1,0.1, 65)) creates a spline table of sin(x) slightly beyond the first period
	"""
	x=[xmin + dx*i for i in xrange(count)]
	return InterpolatingFunction(x,x).SetName('x')

def LogLogInterpolatingGrid(xmin, dx, count):
	"create a log-log interpolating grid with both x & y set to (xmin, xmin*dx, ... xmin * dx**(count -1) )"
	x=[xmin]
	for i in xrange(count-1):
		x.append(x[-1]*dx)
	return LogLogInterpolatingFunction(x,x).SetName('x')


if __name__=="__main__":
	def as(x): return _numeric.array_str(x, precision=3)
	
	ag=ag1=LinearInterpolatingGrid(0.01, 5.0,4)	
	print ag
	
	print C2Constant(11.5)
	print C2Quadratic(x0=5, a=2, b=1, c=0)
	print C2PowerLaw(a=1.5, b=-2.3)
	print LogLogInterpolatingGrid(0.1, 1.1, 20)
	
	print C2Linear(1.3,2.5).apply(ag1)
	print C2Quadratic(0,1,0,0).apply(ag1)
	print ag1*ag1*ag1
		
	fn=C2sin(C2sqrt(ag1*ag1*ag1))	
	print fn
	
	for i in range(7):
		print i, math.sin((i+0.01)**(3./2.)), fn(i+0.01)
	
	x1=fn.find_root(0.1, 1.35128, 0.1, 0.995, trace=True)
	print x1, math.sin(x1**(3./2.)), fn(x1)-0.995
	
	print fn([1., 2., 3., 4.])
	
	
