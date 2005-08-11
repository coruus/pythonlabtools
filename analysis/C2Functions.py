"""A group of classes which make it easy to manipulate smooth functions, including cubic splines. 

C2Functions know how to keep track of the first and second derivatives of functions, and to use this information in, for example, find_root()
to allow much more efficient solutions to problems for which the general solution may be expensive.

The two primary classes are 
	C2Function, which represents an unevaluted function and its derivatives, and 
	InterpolatingFunction,  which represent a cubic spline of provided data.

C2Functions can be combined with unary operators (nested functions) or binary operators (+-*/ etc.)

Developed by Marcus H. Mendenhall, Vanderbilt University Keck Free Electron Laser Center, Nashville, TN USA
email: marcus.h.mendenhall@vanderbilt.edu
Work supported by the US DoD  MFEL program under grant FA9550-04-1-0045
version $Id: C2Functions.py,v 1.28 2005-08-11 21:37:54 mendenhall Exp $
"""
_rcsid="$Id: C2Functions.py,v 1.28 2005-08-11 21:37:54 mendenhall Exp $"

import math
import operator
import types

import Numeric as _numeric #makes it easy to change to NumArray later
_myfuncs=_numeric #can change to math for (possibly greater) speed (?) but no vectorized evaluation

from analysis import spline as _spline #also for later flexibility

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
		
		self.sampling_grid=None
		
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

		ftol=(5e-16*abs(value)+2e-308);
		xtol=(5e-16*(abs(upper_bracket)+abs(lower_bracket))+2e-308);

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
			if disc >= 0:
				if b>=0:
					q=-0.5*(b+math.sqrt(disc))
				else:
					q=-0.5*(b-math.sqrt(disc))

				if q*q > abs(a*c): delta=c/q	#since x1=q/a, x2=c/q, x1/x2=q^2/ac, this picks smaller step
				else: delta=q/a
				root+=delta;

			if disc < 0 or root<lower_bracket or root>upper_bracket:	#if we jump out of the bracket, bisect
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
	def __pow__(self, right):
		"a**b returns a new C2Function which represents the power law, with a optimization for numerical powers"
		return C2Power(self, right)

	def partial_integrals(self, recur_data, **args):
		"""def partial_integrals(self, xgrid, relative_error_tolerance=1e-12, derivs=2, 
			absolute_error_tolerance=1e-12, depth=0, debug=0, extrapolate=1, allow_recursion=True)
			Return the integrals of a function between the sampling points xgrid.  The sum is the definite integral.
			The choices for derivs are 0, 1 or 2, anything else is an error.
			The derivs parameter is used as follows: 
				derivs=0 uses Simpson's rule (no derivative information).   
				derivs=1 uses a 6th order technique based the first derivatives, but no second derivatives
				derivs=2 uses a 9th (really 10th, since the 9th order error vanishes by symmetry) order tehcnique based the first and second derivatives.
				Be very aware that the 9th order method will only really benefit with very smooth functions, but then it is magic!
		"""
	
		if type(recur_data[1]) is not tuple:	#this is an initialization call, set everything up
			funcgrid=[ ( (x, )+ self.value_with_derivatives(x) ) for x in recur_data ] #convert x list to full function list
			self.total_func_evals=len(recur_data)
			
			derivs= args.get('derivs', 2)
			if derivs==0:
				eps_scale, extrap_coef= 0.1, 16.0
			elif derivs==1:
				eps_scale, extrap_coef = 0.1, 64.0
			elif derivs==2:
				eps_scale, extrap_coef = 0.01, 1024.0
			
			else:
				raise C2Exception("Illegal derivs passed to adaptivePartial_integrals, must be 0, 1 or 2, got %d" % derivs)

			allow_rec= args.get('allow_recursion', True)
			extrapolate = args.get('extrapolate', True) and allow_rec #can only extrapolate after recursion
			
			recur_data=[0, funcgrid, ( ),
					args.get('absolute_error_tolerance', 1e-12),
					args.get('relative_error_tolerance', 1e-12),
					args.get('debug', 0),
					extrapolate, derivs, eps_scale, extrap_coef, allow_rec
				]			
		
		depth, funcgrid, old_integrals, absolute_error_tolerance, relative_error_tolerance, debug, extrapolate, derivs, eps_scale, extrap_coef, allow_rec=recur_data
		
		retvals=[0.0]*(len(funcgrid)-1)
		
		for i in range(len(funcgrid)-1):
			x0, y0, yp0, ypp0=funcgrid[i]
			x2, y2, yp2, ypp2=funcgrid[i+1]
			x1=0.5*(x0+x2)
			y1, yp1, ypp1=self.value_with_derivatives(x1)
			self.total_func_evals+=1
			dx=x2-x0
			dx2=0.5*dx

			#check for underflow on step size, which prevents us from achieving specified accuracy. 
			if abs(dx) < abs(x1)*relative_error_tolerance:
				raise C2Exception("Step size underflow in partial_integrals, depth = %d, x = %.4f" % (depth, x1))
			
			#if we are below the top level, the previous term is already computed.  Otherwise, compute it to one lower order
			if  old_integrals:
				total=old_integrals[i]
			elif derivs==2:
				total=( (14*y0 + 32*y1 + 14*y2) +  dx * (yp0 - yp2) ) * dx /60.

			elif derivs==1:
				total=(y0+4.0*y1+y2)*dx/6.0
			else:
				total=0.5*(y0+y2)*dx

			if derivs==2:
				#use ninth-order estimates for each side, from full set of all values (!) (Thanks, Mathematica!)
				left=	( ( (169*ypp0 + 1024*ypp1 - 41*ypp2)*dx2 + (2727*yp0 - 5040*yp1 + 423*yp2) )*dx2 + (17007*y0 + 24576*y1 - 1263*y2) )* (dx2/40320.0)
				right=	( ( (169*ypp2 + 1024*ypp1 - 41*ypp0)*dx2 - (2727*yp2 - 5040*yp1 + 423*yp0) )*dx2 + (17007*y2 + 24576*y1 - 1263*y0) )* (dx2/40320.0)
			elif derivs==1:
				left= 	( (202*y0 + 256*y1 + 22*y2) + dx*(13*yp0 - 40*yp1 - 3*yp2) ) * dx /960.
				right= 	( (202*y2 + 256*y1 + 22*y0) - dx*(13*yp2 - 40*yp1 - 3*yp0) ) * dx /960.
			else:
				left=	(5*y0 + 8*y1 - y2)*dx/24.
				right=	(5*y2 + 8*y1 - y0)*dx/24.

			eps= abs(total-(left+right))*eps_scale #the real error should be 2**order times smaller on this iteration, be conservative
			if extrapolate:
				eps=eps*eps_scale #gain another factor of 2**order if extrapolating (typical), but be conservative
				
			if  (not allow_rec) or eps < absolute_error_tolerance or eps < abs(total)*relative_error_tolerance:
				if not extrapolate:
					retvals[i]=left+right
				else:
					retvals[i]=(extrap_coef*(left+right) - total) / (extrap_coef-1) #since h fell by 2, h**6=64, and we are extrapolating in h**6
				if debug==1: print "accepted results at depth ", depth,  "x, dx = %7.3f, %10.6f" % (x1, dx), "scaled error = ", eps/ (abs(total)*relative_error_tolerance)
			else:
				if debug==1: print "rejected results at depth ", depth,  "x, dx = %7.3f, %10.6f" % (x1, dx), "scaled error = ", eps/ (abs(total)*relative_error_tolerance)
				recur_data[0:4]=depth+1, (funcgrid[i], (x1, y1, yp1, ypp1), funcgrid[i+1]), (left, right), absolute_error_tolerance/2
				l, r =self.partial_integrals(recur_data)
				retvals[i]=l+r
		
		if debug ==2:
			print "\n integrating at depth ", depth
			print xgrid
			print funcgrid
			import operator
			print retvals
			if old_integrals: 
				print map(operator.sub, old_integrals, retvals)
			print "\n returning from depth ", depth
		
		return retvals
	
	def SetSamplingGrid(self, grid):
		"establish a set of 'interesting' points for starting to sample this function"
		self.sampling_grid=_numeric.array(grid)
	
	def GetSamplingGrid(self, xmin, xmax):
		"get a set of points, including xmin and xmax, which are reasonable points to evaluate the function between them"
		if not self.sampling_grid or xmin > self.sampling_grid[-1] or xmax < self.sampling_grid[0]:
			return (xmin, xmax) #don't really have any information on the funciton in this interval.
		else:
			sg=self.sampling_grid
			
			firstindex, lastindex=_numeric.searchsorted(sg, (xmin, xmax))
			lastindex=lastindex-1
						
			grid=[xmin]+list(sg[firstindex:lastindex+1])+[xmax] #insert points  from source grid to destination into the middle of our grid
			
			if len(grid) > 2 : #attempt to clear out points put too close together at the beginning, if we have at least three points
				x0, x1, x2 = grid[:3]
				ftol=10.0*(1e-14*(abs(x0)+abs(x1))+1e-300)
				if (x1-x0) < ftol or (x1-x0) < 0.1*(x2-x0): del grid[1] #remove collision
			
			if len(grid) > 2 : #attempt to clear out points put too close together at the end, if we have at least three points
				x0, x1, x2 = grid[-3:]
				ftol=10.0*(1e-14*(abs(x0)+abs(x2))+1e-300)
				if (x2-x1) < ftol or (x2-x1) < 0.1*(x2-x0): del grid[-2] #remove collision

			return grid

	def integral(self, xmin, xmax, **args):
		"carry out integration using our sampling grid"
		grid=self.GetSamplingGrid(xmin, xmax)
		return sum(self.partial_integrals(grid, **args))

	def NormalizedFunction(self, xmin, xmax, norm=1):
		"return a function whose integral on [xmin, xmax] is norm"
		intg=self.integral(xmin, xmax)
		return C2ScaledFunction(self, norm/intg)

	def SquareNormalizedFunction(self, xmin, xmax, weight=None, norm=1):
		"return a function whose square integral on [xmin, xmax] with optional weight function is norm"
		mesquared=C2Quadratic(a=1.0)(self)
		if weight is not None:
			mesquared=mesquared*weight
			
		grid=self.GetSamplingGrid(xmin, xmax)
		intg=sum(mesquared.partial_integrals(grid))
		return C2ScaledFunction(self, math.sqrt(norm/intg))
	
class C2ScaledFunction(C2Function):
	"a lightweight class to create a function scaled vertically by a scale factor"
	ClassName='Scaled'
	def __init__(self, fn, yscale):
		C2Function.__init__(self)
		self.fn=fn
		self.yscale=yscale
		self.name=fn.name+'* %g' % yscale
	
	def value_with_derivatives(self, x): 
		y, yp, ypp = fn.value_with_derivatives(x)
		ys=self.yscale
		return y*ys, yp*ys, ypp*ys

class C2Constant(C2Function):
	"a constant and its derivatives"
	ClassName='Constant'
	def __init__(self, val):
		C2Function.__init__(self)
		self.val=val
		self.name='%g' % val

	def reset(self, val):
		"reset the value to a new value"
		self.val=val

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
	
	def reset(self, slope=None, y0=None):
		"reset the slope or intercept to a new value"
		if slope is not None:
			self.slope=slope
		if y0 is not None:
			self.y0=y0
			
	def value_with_derivatives(self, x):
		return x*self.slope+self.y0, self.slope, 0.0

class C2Quadratic(C2Function):
	"a*(x-x0)**2 + b*(x-x0) + c"
	def __init__(self, x0=0.0, a=1.0, b=0.0, c=0.0):
		C2Function.__init__(self)
		self.x0, self.a, self.b, self.c = x0, a, b, c
		self.name="(%g*(x-x0)^2 + %g*(x-x0) + %g, x0=%g)" % (a, b, c, x0)
		
	def reset(self, x0=None, a=None, b=None, c=None):
		"reset the parameters to new values"
		if x0 is not None:
			self.x0=x0
		if a is not None:
			self.a=a
		if b is not None:
			self.b=b
		if b is not None:
			self.c=c

	def value_with_derivatives(self, x):
		dx=x-self.x0
		return self.a*dx*dx+self.b*dx+self.c, 2*self.a*dx+self.b, 2*self.a

class C2PowerLaw(C2Function):
	"a*x**b where a and b are constant"
	def __init__(self, a=1.0, b=2.0):
		C2Function.__init__(self)
		self.a, self.b = a, b
		self.b2=b-2
		self.bb1=b*(b-1)
		self.name='%g*x^%g' % (a,b)
		
	def reset(self, a=None, b=None):
		"reset the parameters to new values"
		if a is not None:
			self.a=a
		if b is not None:
			self.b=b

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
		self.left=self.convert_arg(left)
		self.right=self.convert_arg(right)
		C2Function.__init__(self, self.left, self.right)
		
		if isinstance(left, C2BinaryFunction): p1, p2 = '(', ')'
		else: p1, p2='', ''
		self.name=p1+self.left.name+p2+self.name+self.right.name #put on parentheses to kepp hierachy obvious
		
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

class C2Power(C2BinaryFunction):
	"C2power(a,b) returns a new C2Function which evaluates as a^b where neither is necessarily constant.  Checks if b is constant, and optimizes if so"
	name='^'
	def combine(self, x):
		y0, yp0, ypp0=self.left.value_with_derivatives(x)
		y1, yp1, ypp1=self.right.value_with_derivatives(x)
		if isinstance(self.right, C2Constant): #all derivatives of right side are zero, save some time
			ab2=_myfuncs.power(y0, (y1-2.0)) #this if a^(b-2), which appears often
			ab1=ab2*y0 #this is a^(b-1)
			ab=ab1*y0 #this is a^b
			ab1 *= yp0
			ab1 *= y1 #ab1 is now the derivative
			ab2 *= y1*(y1-1)*yp0*yp0+y0*y1*ypp0 #ab2 is now the second derivative
			return ab,  ab1, ab2
		else:
			loga=_myfuncs.log(y0)
			ab2=_myfuncs.exp(loga*(y1-2.0)) #this if a^(b-2), which appears often
			ab1=ab2*y0 #this is a^(b-1)
			ab=ab1*y0 #this is a^b
			yp=ab1*(yp0*y1+y0*yp1*loga)
			ypp=ab2*(y1*(y1-1)*yp0*yp0+2*y0*yp0*yp1*(1.0+y1*loga)+y0*(y1*ypp0+y0*(ypp1+yp1*yp1*loga)*loga))
			return ab, yp, ypp
	
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
		self.xInverted=False
		
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
			self.fXin, self.fXinP, self.fXinPP, self.fXout = self.XConversions =  lambda x: x, lambda x: 1, lambda x: 0, lambda x: x
			self.xNonLin=False
			self.X=_numeric.array(x)
		else:
			self.fXin, self.fXinP, self.fXinPP, self.fXout = self.XConversions
			self.xNonLin=True
			self.X=_numeric.array([self.fXin(q) for q in x])
			if lowerSlope is not None: lowerSlope /= self.fXinP(x[0])
			if upperSlope is not None: upperSlope /= self.fXinP(x[-1])

			if self.X[0] > self.X[-1]: #make sure our transformed X array is increasing
				self.Xraw=self.Xraw[::-1]
				self.X=self.X[::-1]
				self.F=self.F[::-1]
				self.xInverted=True
				lowerSlope, upperSlope=upperSlope, lowerSlope
			
		if not self.xInverted:
			self.SetLowerExtrapolation=self.SetLeftExtrapolation
			self.SetUpperExtrapolation=self.SetRightExtrapolation
		else:
			self.SetLowerExtrapolation=self.SetRightExtrapolation
			self.SetUpperExtrapolation=self.SetLeftExtrapolation
		
		dx=self.X[1:]-self.X[:-1]
		if min(dx) <  0 or min(self.X) < self.X[0] or max(self.X) > self.X[-1]:
			raise C2Exception("monotonicity error in X values for interpolating function: "  + 
				_numeric.array_str(self.X))
			

		self.y2=_spline.spline(self.X, self.F, yp1=lowerSlope, ypn=upperSlope)

	def value_with_derivatives(self, x): 
		if self.xNonLin or self.yNonLin: #skip this if this is a completely untransformed spline
			y0, yp0, ypp0=_spline.splint(self.X, self.F, self.y2, self.fXin(x), derivs=True)
			y=self.fYout(y0)
			fpi=1.0/self.fYinP(y)
			gp=self.fXinP(x)

			# from Mathematica Dt[InverseFunction[f][g[y[x]]], x]
			yprime=gp*yp0*fpi # the real derivative of the inverse transformed output 
			fpp=self.fYinPP(y)
			gpp=self.fXinPP(x)
			#also from Mathematica Dt[InverseFunction[f][g[y[x]]], {x,2}]
			yprimeprime=(gp*gp*ypp0 + yp0*gpp - gp*gp*yp0*yp0*fpp*fpi*fpi)*fpi; 
			return y, yprime, yprimeprime
		else:
			return _spline.splint(self.X, self.F, self.y2, x, derivs=True)

	def SetLeftExtrapolation(self, bound):
		"""Set extrapolation on left end of data set.  
		This will be dynamically assigned to either SetLowerExtrapolation or SetUpperExtrapolation by the constructor
		"""
		xmin=self.fXin(bound)
		if xmin >= self.X[0]: 
			raise C2Exception("Attempting to extrapolate spline within its current range: bound = %g, bounds [%g, %g]" % 
				((bound,) + self.GetDomain()) )
		
		self.X, self.F, self.y2=_spline.spline_extension(self.X, self.F, self.y2, xmin=xmin)
		self.Xraw=_numeric.concatenate(((bound, ), self.Xraw))
		self.SetDomain(min(self.Xraw), max(self.Xraw))
		
	def SetRightExtrapolation(self, bound):
		"""Set extrapolation on right end of data set.  
		This will be dynamically assigned to either SetLowerExtrapolation or SetUpperExtrapolation by the constructor
		"""
		xmax=self.fXin(bound)
		if xmax <= self.X[-1]: 
			raise C2Exception("Attempting to extrapolate spline within its current range: bound = %g, bounds [%g, %g]" % 
				((bound,) + self.GetDomain()) )

		self.X, self.F, self.y2=_spline.spline_extension(self.X, self.F, self.y2, xmax=xmax)
		self.Xraw=_numeric.concatenate((self.Xraw, (bound, )))
		self.SetDomain(min(self.Xraw), max(self.Xraw))

	def YtoX(self):
		"returns a new InterpolatingFunction with our current grid of Y values as the X values"
		
		yv=self.fYout(self.F) #get current Y values transformed out
		if yv[1] < yv[0]: yv=yv[::-1]
		f=InterpolatingFunction(yv, yv, XConversions=self.XConversions, YConversions=self.YConversions)
		f.SetName("x values: "+self.name)
		return f
		
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

LogConversions=_myfuncs.log, lambda x: 1.0/x, lambda x: -1.0/(x*x), _myfuncs.exp

class LogLinInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores log(x) vs. y"
	ClassName='LogLinInterpolatingFunction'
	XConversions=LogConversions

class LinLogInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores x vs. log(y), useful for functions with exponential-like behavior"
	ClassName='LinLogInterpolatingFunction'
	YConversions=LogConversions

class LogLogInterpolatingFunction(InterpolatingFunction):
	"An InterpolatingFunction which stores log(x) vs. log(y), useful for functions with power-law-like behavior"
	ClassName='LogLogInterpolatingFunction'
	XConversions=LogConversions
	YConversions=LogConversions
	def log_log_partial_integrals_broken(self, xgrid, debug=False):
		"""Return the integrals of a function between the sampling points xgrid.  The sum is the definite integral.
		This version knows about integrating in log-log space, at least to some simple approximation.
		The approximation is that ln(y)=a + b*ln(x) + c*ln(x)^2 -> a + (b + c*ln(xbar)) * ln(x) where xbar=sqrt(x0*x1) for the interval 
		"""
		
		xgridl=_numeric.log(xgrid)
		dx=xgridl[1:]-xgridl[:-1]
	
		#compute all log values & derivatives of logs at sampling points
		y, yp, ypp=_spline.splint(self.X, self.F, self.y2, xgridl, derivs=True) 
		
		localpower=yp[:-1]+0.25*ypp[:-1]*(xgridl[1:]+xgridl[:-1])+1 #this is the 1+slope of log x assuming  ln^2(x) is ln(x)*ln(sqrt(x0*x1))
		if debug: 
			print "\ndebug log_log_partial_integrals:"
			print dx
			print y
			print yp
			print ypp
			print localpower
			print zip(self.fYout(y[:-1]), self.fYout(localpower*dx))
		partials=_numeric.exp(y[:-1]-(localpower+1)*xgridl[:-1])*(_numeric.exp(localpower*dx)-1.0)/localpower
		
		return partials

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

class AccumulatedHistogram(InterpolatingFunction):
	"""An InterpolatingFunction which is the cumulative integral of the (histogram) specified by binedges and binheights.
		Note than binedges should be one element longer than binheights, since the lower & upper edges are specified. 
		Note that this is a malformed spline, since the second derivatives are all zero, so it has less continuity.
		Also, note that the bin edges can be given in backwards order to generate the reversed accumulation (starting at the high end) 
	"""
	ClassName='AccumulatedHistogram'

	def __init__(self, binedges, binheights, normalize=False, inverse_function=False, drop_zeros=True, **args):
		be=_numeric.array(binedges, _numeric.Float)
		bh=_numeric.array(binheights, _numeric.Float)

		if drop_zeros or inverse_function: #invese functions cannot have any zero bins or they have vertical sections
		
			nz=_numeric.not_equal(bh, 0) #mask of non-empty bins, or lower edges
			if not inverse_function: nz[0]=nz[-1]=1 #always preserve end bins to keep X range, but don't dare do it for inverses
			bh=_numeric.compress(nz, bh)
			be=_numeric.compress(_numeric.concatenate( (nz, (nz[-1],) ) ), be)
				
		cum=_numeric.concatenate( ( (0,), _numeric.cumsum( (be[1:]-be[:-1])*bh ) ))
		
		if be[1] < be[0]: #fix backwards bins, if needed.
			be=be[::-1]
			cum=cum[::-1]
			cum*=-1 #the dx values were all negative if the bins were backwards, so fix the sums

		if normalize:
			cum *= (1.0/max(cum[0], cum[-1])) #always normalize on the big end

		if inverse_function:
			be, cum = cum, be #build it the other way around
			
		InterpolatingFunction.__init__(self, be, cum, **args)
		self.y2 *=0 #clear second derivatives... we know nothing about them

class LogLogAccumulatedHistogram(AccumulatedHistogram):
	"same as AccumulatedHistogram, but log-log axes"
	ClassName='LogLogAccumulatedHistogram'
	XConversions=LogConversions
	YConversions=LogConversions

class InverseIntegratedDensity(InterpolatingFunction):
	"""InverseIntegratedDensity starts with a probability density function, generates the integral, 
	and generates a LinLogInterpolatingFunction which, when evaluated using a uniform random on [0,1] returns values
	with a density distribution equal to the input distribution
	If the data are passed in reverse order (large X first), the integral is carried out from the big end,
	and then the data are reversed to the result in in increasing X order.  
	"""
	
	IntermediateInterpolator=InterpolatingFunction
	
	def __init__(self, bincenters, binheights):
		
		if not isinstance(binheights, C2Function): #must be a data table, create interpolator
			np=len(binheights)	
			be=_numeric.array(bincenters)
			bh=_numeric.array(binheights)
				
			reversed = be[1] < be[0]	#// check for backwards  channels
		
			if reversed:
				be=be[::-1]
				bh=bh[::-1]
		
			temp=self.IntermediateInterpolator(be, bh)		#// create a temporary InterpolatingFunction to integrate
		else:
			temp=binheights
			
		# integrate from first to last bin in original order, leaving results in integral
		# ask for relative error of 1e-6 on each bin, with absolute error set to 0 (since we don't know the data scale).
		integral=[0] + temp.partial_integrals(bincenters, 
				absolute_error_tolerance=0, 
				relative_error_tolerance=1e-6) 
	
		scale=1.0/sum(integral)	

		for i in range(1,len(integral)):
			integral[i]=integral[i]*scale + integral[i-1]
		integral[-1]=1.0 #force exact value on the boundary	
		
		InterpolatingFunction.__init__(self, integral, bincenters, 
							lowerSlope=1.0/(scale*temp(bincenters[0])), upperSlope=1.0/(scale*temp(bincenters[-1]))
				) 	# use integral as x axis in inverse function

class LinLogInverseIntegratedDensity(InverseIntegratedDensity):
	YConversions=LogConversions
	IntermediateInterpolator=LogLogInterpolatingFunction

if __name__=="__main__":
	print _rcsid
	def as(x): return _numeric.array_str(x, precision=3)
	
	if 0:
		ag=ag1=LinearInterpolatingGrid(1, 1.0,4)	
		print ag	
		try:
			ag.SetLowerExtrapolation(2)
		except:
			import sys
			print "***got expected error on bad extrapolation: ", sys.exc_value
		else:
			print "***failed to get expected exception on bad extrapolation"
			
		ag.SetLowerExtrapolation(-2)
		ag.SetUpperExtrapolation(15)
		print ag
		
		print C2Constant(11.5)
		print C2Quadratic(x0=5, a=2, b=1, c=0)
		print C2PowerLaw(a=1.5, b=-2.3)
		print LogLogInterpolatingGrid(0.1, 1.1, 20)
		
		print C2Linear(1.3,2.5).apply(ag1)
		print C2Quadratic(0,1,0,0).apply(ag1)
		print ag1*ag1*ag1
		print (ag1*ag1*ag1).YtoX()
			
		try:
			ag13=(ag1*ag1).YtoX()
		except:
			import sys
			print "***got expected error on bad X axis: ", sys.exc_value
		else:
			print "***failed to get expected exception on bad X axis"
			
		fn=C2sin(C2sqrt(ag1*ag1*ag1)).SetDomain(0, ag1.GetDomain()[1]) #cut off sqrt(negative)
		print fn
		
		for i in range(10):
			print i, "%20.15f %20.15f" % (math.sin((i+0.01)**(3./2.)), fn(i+0.01) )
		
		x1=fn.find_root(0.0, 1.35128, 0.1, 0.995, trace=True)
		print x1, math.sin(x1**(3./2.)), fn(x1)-0.995
		
		print fn([1., 2., 3., 4.])
		
		print "\nPower law sin(x)**log(x) tests"
		fn=C2sin**C2log
		for i in range(10):
			x=0.1*i + 6.4
			print ( "%10.3f "+3*"%15.12f ")%( (x, )+fn.value_with_derivatives(x))
			
		print "\nPower law sin(x)**3.2 tests"
		fn=C2sin**3.2
		for i in range(10):
			x=0.1*i + 6.4
			print ( "%10.3f "+3*"%15.12f ")%( (x, )+fn.value_with_derivatives(x))

		import math
		print "\nIntegration tests"
		sna=C2sin(C2PowerLaw(1,2)) #integrate sin(x**2)
		for sample in (5, 11, 21, 41, 101):
			xg=_numeric.array(range(sample), _numeric.Float)*(2.0/(sample-1))
			partials=sna.partial_integrals(xg)
			if sample==10: print _numeric.array_str(partials, precision=8, suppress_small=False, max_line_width=10000)
			sumsum=sum(partials)
			yg=sna(xg)
			simp=sum(sna.simpson_partial_integrals(xg[::2])) #sparsify grid by 2 to be fair since Simpsons fills in extra points
			exact=0.804776489343756110
			print ("%3d "+6*"%18.10g") %  (sample, 
				simp, simp-exact, (exact-simp)*sample**4, 
				sumsum, sumsum-exact, (exact-sumsum)*sample**6) #the comparision is to the Fresnel Integral from Mathematica
	
	if 0:
		print "\nBessel Functions by integration"
		print """Warning... the adaptive integrator looks worse than the non-adaptive one here. 
		There is some subtle cancellation which makes uniform sampling give extremely accurate answers for Bessel's integral, 
		so it isn't the fault of the adaptive integrator.  The simple one works way too well, here, by accident.
		"""
			
		def bessj(n, z, point_density=2):
			f=C2cos(C2Linear(slope=n) - C2Constant(z) * C2sin)
			pc=int((abs(z)+abs(n)+2)*point_density)
			g=_numeric.array(range(pc), _numeric.Float)*(math.pi/(pc-1))
			return sum(f.partial_integrals(g, allow_recursion=False))/math.pi, f.total_func_evals
		
		def bessj_adaptive(n, z, derivs=1):
			f=C2cos(C2Linear(slope=n) - C2Constant(z) * C2sin)
			pc=8
			g=_numeric.array(range(pc), _numeric.Float)*(math.pi/(pc-1))
			return sum(f.partial_integrals(g, absolute_error_tolerance=1e-14, relative_error_tolerance=1e-14, debug=0, derivs=derivs))/math.pi, f.total_func_evals

		for order, x in ( (0, 0.1), (0,5.5), (2,3.7), (2,30.5)):
			v1, n1=bessj(order, x)
			v2, n2=bessj_adaptive(order, x, 1)
			v3, n3=bessj_adaptive(order, x, 2)
			print ("%6.2f %6.2f %6d %6d %6d "+3*"%20.15f ") % (order, x, n1, n2, n3, v1, v2, v3 )
		

	if 0:
		print "\nLogarithms  by integration"
		
		pc=3
		for lv in (0.1, 1.0, 2.0, 5.0, 10.0):
			b=math.exp(lv)
			np=int(pc*b)+4
			g=_numeric.array(range(np), _numeric.Float)*(b-1)/(np-1)+1

			v0=C2recip.partial_integrals(g, allow_recursion=False)
			n0=C2recip.total_func_evals
			
			
			v1=C2recip.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)), absolute_error_tolerance=1e-6, 
					extrapolate=1, debug=0, derivs=0)
			n1=C2recip.total_func_evals

			v2=C2recip.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)), absolute_error_tolerance=1e-12, 
					extrapolate=1, debug=0, derivs=1)
			n2=C2recip.total_func_evals
			
			v3=C2recip.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)), absolute_error_tolerance=1e-12, 
					extrapolate=1, debug=0, derivs=2)
			n3=C2recip.total_func_evals

			print ("%20.15f %10.2f %6d %6d %6d %6d "+4*"%20.15f ") % (lv, b, n0, n1, n2, n3, sum(v0), sum(v1), sum(v2), sum(v3) )

			
	if 0:
		print "\nPowers  by integration"
		
		fn=C2recip(C2Quadratic(a=1.0, b=0.01, c=0.01)) #make approximate power law
		
		pc=5
		for lv in (0.1, 1.0, 2.0, 5.0, 10.0):
			b=1.0+lv
			np=int(pc*b)+4
			g=_numeric.array(range(np), _numeric.Float)*(b-1)/(np-1)+1

			v0=fn.partial_integrals(g, allow_recursion=False)
			n0=fn.total_func_evals
						
			v1=fn.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)),
					absolute_error_tolerance=1e-8,  relative_error_tolerance=1e-8,
					extrapolate=0, debug=0, derivs=0)
			n1=fn.total_func_evals

			v2=fn.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)),
					absolute_error_tolerance=1e-14,  relative_error_tolerance=1e-14,
					extrapolate=1, debug=0, derivs=1)
			n2=fn.total_func_evals
		
			v3=fn.partial_integrals(_numeric.array((1.0,math.sqrt(b), b)),
					absolute_error_tolerance=1e-14,  relative_error_tolerance=1e-14,
					extrapolate=1, debug=0, derivs=2)
			n3=fn.total_func_evals

			print ("%20.15f %10.2f %6d %6d %6d %6d "+4*"%20.15f ") % (lv, b, n0, n1, n2, n3, sum(v0), sum(v1), sum(v2), sum(v3) )

	if 0:
		print "\nAccumulatedHistogram tests"
		xg=(_numeric.array(range(21), _numeric.Float)-10.0)*0.25
		yy=_numeric.exp(-xg[:-1]*xg[:-1])
		yy[3]=yy[8]=yy[9]=0
		ah=AccumulatedHistogram(xg[::-1], yy[::-1], normalize=True)
		print ah([-2, -1, 0, 1, 2])
		ah=AccumulatedHistogram(xg[::-1], yy[::-1], normalize=True, drop_zeros=False)
		print ah([-2, -1, 0, 1, 2])
		ahi=AccumulatedHistogram(xg, yy, normalize=True, inverse_function=True)
		print ahi([0, 0.01,0.5, 0.7, 0.95, 1])
	
	
	if 0:
		xv=_numeric.array([1.5**(0.1*i) for i in range(100)])
		yv=_numeric.array([x**(-4)+0.25*x**(-3) for x in xv])
		f=LogLogInterpolatingFunction(xv,yv, 
			lowerSlope=-4*xv[0]**(-5)+(0.25)*(-3)*xv[0]**(-4),
			upperSlope=-4*xv[-1]**(-5)+(0.25)*(-3)*xv[-1]**(-4))
		f0=C2PowerLaw(1., -4)+C2PowerLaw(0.25, -3)
		print f(_numeric.array([1.,2.,3.,4., 5., 10., 20.]))
		print f0(_numeric.array([1.,2.,3.,4., 5., 10., 20.]))
		partials=f.partial_integrals(_numeric.array([1.,2.,3.,4., 5., 10., 20., 21.]), debug=False)
		print partials, sum(partials)
		partials=f.log_log_partial_integrals(_numeric.array([1.,2.,3.,4., 5., 10., 20., 21.]), debug=False)
		print partials, sum(partials)
		pp=f0.partial_integrals(_numeric.array(range(11), _numeric.Float)*0.1 + 20)
		print pp, sum(pp)
	
	if 0:
		print "\nTesting LinLogInverseIntegratedDensity for 1/e^2"
		energies=[float(2**(0.5*n)) for n in range(41)]
		spect=[10000.0/(e*e) for e in energies]
		
		e0=energies[-1]
		e1=energies[0]
		
		pf=LinLogInverseIntegratedDensity(energies[::-1], spect[::-1])
		
		print energies[0], energies[-1]
		
		for i in range(41):
			r=(0.025*i)**2
			mma=(e0*e1)/(r*e0 + e1 - r*e1)
			print "%12.6f %12.4e %12.4e %12.4f " % (r, pf(r), mma, 100*(pf(r)-mma)/mma)
		
		print "\nTesting LinLogInverseIntegratedDensity for 1/e using a C2Function instead of a table"
		energies=[float(2.0**n) for n in range(21)]
		
		e0=energies[-1]
		e1=energies[0]
				
		pf=LinLogInverseIntegratedDensity(energies[::-1], C2PowerLaw(a=1000.0, b=-1))
		
		print energies[0], energies[-1]
		
		for i in range(41):
			r=(0.025*i)**2
			mma=e0*(e1/e0)**r
			
			print "%12.6f %12.4e %12.4e %12.4e " % (r, pf(r), mma, 100*(pf(r)-mma)/mma)

	if 1:
		fn=C2sin *2.0 #make a new function
		print "\nInitial sampling grid =", 
		grid=(0., 3., 6., 9., 12.)
		for i in range(len(grid)): print grid[i],
		print
		
		fn.SetSamplingGrid(grid)
		
		print "Starting tests of samples from grid"
		
		for xmin,xmax in ( (-10,-1), (20,30), (-3, 15), (-3, 2), (-3, 7), (5.0, 5.5), (1., 7.), (5.99, 10), (2, 9.01), (5., 20.) ):
			v=fn.GetSamplingGrid(xmin,xmax)
			print "%10.3f %10.3f: " %(xmin, xmax),  
			for i in range(len(v)): print v[i],
			print
			
		sn=fn.NormalizedFunction(0., math.pi)
		
		print "integral of non-normalized and normalized function ", fn.integral(0., math.pi), sn.integral(0., math.pi)
		
		gn=fn.SquareNormalizedFunction(0., 4.0*math.pi)
		
		fn2=fn*fn
		fn2.SetSamplingGrid((0., math.pi, 2.0*math.pi, 3.0*math.pi, 4.0*math.pi))
		gn2=gn*gn
		gn2.SetSamplingGrid((0., math.pi, 2.0*math.pi, 3.0*math.pi, 4.0*math.pi))
		
		print "integral of square of non-normalized and square-normalized function ", fn2.integral(0., 4.0*math.pi), gn2.integral(0., 4.0*math.pi)
		
	
