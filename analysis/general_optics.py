from math import *
import math
import Numeric
import types
import traceback
import LinearAlgebra

clight=299792458. #m/sec
deg=pi/180.0

Infinity="Infinity"

ambient_index=1.0003 #air

def get_ambient_index(lambda0):
	if type(ambient_index) is types.FunctionType:
		return ambient_index(lambda0)
	else:
		return ambient_index

def reset():
	global ambient_index
	ambient_index=1.0003

def vec_mag(vec):
	return math.sqrt(Numeric.sum(vec*vec))

def planesolve(x0, u, x1, v):
	"planesolve(x0, u, x1, v) returns a, x such that x_=a*u_+x0_ solves (x_-x1_).v_=0"
	if vec_mag(x1-x0) < 1e-10:
		return 0, x1
	a=Numeric.dot((x1-x0), v) / Numeric.dot(u,v)
	#print "\nplanesolve: ", x0, u, x1, v, a, a*u+x0, "\n"
	return a, a*u+x0

def cross(a,b):
	x1,y1,z1=tuple(a)
	x2,y2,z2=tuple(b)
	return Numeric.array((y1*z2-z1*y2, z1*x2-x1*z2, x1*y2-y1*x2))
	
def wrap_angle(theta):
	"wrap_angle always returns an angle in (-180,180)"
	return ( (theta + 180) % 360 ) -180

def sincosdeg(theta_deg):
	return math.sin(theta_deg*deg), math.cos(theta_deg*deg)

def normalize(x,y):
	rho=math.sqrt(x*x+y*y)
	if rho < 1e-10:
		return 1., 0.
	else:
		return x/rho, y/rho

def eulermat(theta, eta, phi):
	"eulermat(theta, eta, phi) returns r(eta)*r(theta)*r(phi) phi=orientation around normal, theta=yaw, eta=normal roll"
	pm=Numeric.identity(3).astype(Numeric.Float)
	em=Numeric.identity(3).astype(Numeric.Float)
	tm=Numeric.identity(3).astype(Numeric.Float)
	sp, cp=sincosdeg(phi)
	st, ct=sincosdeg(theta)
	se, ce=sincosdeg(eta)
	pm[0,0]=pm[1,1]=cp; pm[0,1]=-sp; pm[1,0]=sp
	em[0,0]=em[1,1]=ce; em[0,1]=-se; em[1,0]=se
	tm[0,0]=tm[2,2]=ct; tm[2,0]=-st; tm[0,2]=st
	return Numeric.dot(em, Numeric.dot(tm, pm))


class euler:
	def __init__(self, theta, eta, phi):
		self.theta=theta
		self.eta=eta
		self.phi=phi
		self.mat=eulermat(theta, eta, phi)
	
	def __str__(self):
		s=""
		if self.theta != 0:
			s=s+"theta = %.1f, " % self.theta
		if self.eta != 0:
			s=s+"eta = %.1f, " % self.eta
		if self.phi != 0:
			s=s+"phi = %.1f" % self.phi
		if s[-2:]==", ": s=s[:-2]
		return s

def general_simil_tens(transform, tens, reversed=0):
	"return transform^(-1)*tens*transform or the other way around if reversed=1"
	q,r,s,t=transform[0,0], transform[0,1], transform[1,0], transform[1,1]
	inverse=Numeric.array(((t,-r),(-s,q))) / (q*t-r*s)
	if reversed:
		inverse, transform=transform, inverse
	return Numeric.dot(inverse, Numeric.dot(tens, transform))

def simil_tens_cs(k, s, tens):
	"similarity transform a 2x2 tensor by angle k=cos(theta), s=sin(theta)."
	a,b,c,d=tens[0,0], tens[0,1], tens[1,0], tens[1,1]
	s2=s*s; k2=k*k; sk=s*k
	a1=a*k2+d*s2+sk*(b+c)
	b1=b*k2-c*s2+sk*(d-a)
	c1=c*k2-b*s2+sk*(d-a)
	d1=d*k2+a*s2-sk*(b+c)
	return Numeric.array(((a1,b1),(c1,d1)))

def simil_tens(theta, tens):
	"similarity transform a 2x2 tensor.  If the tensor is deeper than 2nd rank, the inner indices are passed through"
	return simil_tens_cs(math.cos(theta), math.sin(theta), tens)

def principal_axis_angle(q22):
	a,b,c,d=tuple(Numeric.ravel(q22))
	if abs(d-a) < 1e-10: return 0
	#if we are passed a q-tensor for which 2b/(d-a) isn't real, this is wrong!
	y=(2*b/(d-a))
	qtheta=math.atan(y.real)/2.0 #principal axis direction of q
	if abs(y.real) > 1e-20 and abs(y.imag)/abs(y.real) > 1e-6:
		raise "Bad tensor to diagonalize: y="+str(y)+" theta=: "+str(qtheta/deg)+"\n"+str(q22)
		
	return qtheta

def expand_to_2x2tensor(object):
	a=Numeric.array(object)
	if a.shape==():
		return a[0]*Numeric.identity(2)
	elif a.shape==(2,):
		return Numeric.array(((a[0],0),(0,a[1])))
	else:
		return a
		

class qtens:
	def __init__(self, lambda_initial, q=None, w=None, r=None, name=None, qit=None, medium_index=1.0):
		self.name=name			
		self.lambda0=lambda_initial*medium_index
		self.medium_index=medium_index
		if qit is not None: #we are given a real inverse-q tensor, just use it
			self.qit=qit
		else:
			if q is None:
				if r is Infinity:
					rinv=0
				else:
					rinv=1.0/r
				qi=complex(rinv , -self.lambda0 / (self.medium_index*(math.pi*w**2)) )
				self.qit=Numeric.array(((qi,0),(0,qi)))
			elif type(q)==types.ComplexType:
				self.qit=Numeric.array(((1.0/q,0),(0,1.0/q)))
			else:
				self.qit=Numeric.array(((1.0/q[0],0),(0,1.0/q[1])))
				
	def q_moments(self):
		u,v=LinearAlgebra.eigenvectors(self.qit)
		return v, 1.0/u[0], 1.0/u[1]
	
	def qi_moments(self):
		u,v = LinearAlgebra.eigenvectors(self.qit)
		return v, u[0], u[1]


	def set_medium_index(self, new_index):
		self.qit=self.qit*(self.medium_index/new_index)
		self.medium_index=new_index
		
	def rw(self, qi):
		w=math.sqrt(-self.lambda0/qi.imag/math.pi/self.medium_index)
		if abs(qi.real) < 1e-15:
			r=Infinity
		else: 
			r=1.0/qi.real			
		return r, w
	
	def next_waist(self):
		theta, qx, qy=self.q_moments()
		qxr=qx.real
		qyr=qy.real
		if qxr*qyr < 0:
			return -min((qxr,qyr)) #one is negative, use it
		else:
			return -max((qxr, qyr)) #the least negative one is the first waist, or no waist if both positive
		
	def rwstr(self, qi):
		r,w=self.rw(qi)
		if r is not Infinity:
			r="%.2f" % r
		dz=(1.0/qi).real
		return ("r=%s w=%.3e dz=%.3f" % (r,w, dz)) # +" qi=("+str(qi)+")"
		
	def __str__(self):
		theta, q0, q1=self.qi_moments()
		if self.name is not None:
			s=self.name+": "
		else: 
			s=""
			
		if abs(q0-q1) < 1e-8:
			s=s+"q={"+self.rwstr(q0)+"}"
		else:
			s=s+"qxx={"+self.rwstr(q0)+"}, qyy={"+self.rwstr(q1)+"}"
			
		return s
	def __repr__(self):
		return self.__str__()
		
	def qw(self, element):
		element.q_transform(self)		
		return self #make daisy-chaining easy
	
	def abcd_transform(self, abcd):
		dot=Numeric.dot
		ar=Numeric.array
		a, b, c, d = [expand_to_2x2tensor(i) for i in ar(abcd).flat]
		q,r,s,t=(a+dot(b, self.qit)).flat
		inverse=ar(((t,-r),(-s,q))) / (q*t-r*s)
		self.qit=dot(c+dot(d, self.qit), inverse)
	
	def drift(self, length):
		q,r,s,t=(Numeric.identity(2)+length*self.qit).flat
		inverse=Numeric.array(((t,-r),(-s,q))) *(1.0/ (q*t-r*s))
		self.qit=Numeric.dot(self.qit, inverse)

	def focus(self, strength):
		self.qit=self.qit+expand_to_2x2tensor(strength)

	def clone(self, newname=None):
		q=copy.deepcopy(self)
		if newname is not None:
			q.name=newname
		return q

	def transform(self, tensor):
		tr=Numeric.transpose
		dot=Numeric.dot
		self.qit=dot(tensor, dot(self.qit, Numeric.conjugate(tr(tensor))))
		return self
		
import exceptions

class OpticDirectionError(exceptions.AssertionError):
	"Light pointing the wrong way to hit this optic!"
import copy

class beam:
	"the beam class carries around information for a beam of light: starting point, polarization, wavelength... starting along +z"
	def __init__(self, x0, q, lam=None, polarization=(1,0), medium_index=None):
		
		self.x0=Numeric.array(x0,Numeric.Float)
		self.matrix_to_global=Numeric.identity(3).astype(Numeric.Float)
		self.polarization=Numeric.array(polarization, Numeric.Complex)
		if isinstance(q, qtens):
			self.q=q.clone()
		else:
			if medium_index is None:
				medium_index=get_ambient_index(lam)
			self.q=qtens(lam, q=q, medium_index=medium_index)
		self.total_drift=0.0 #accumulated travel distance
		self.total_drift_time=0.0 #accumulated distance*index
		self.marks={"order":[]}
	
	def direction(self):
		return Numeric.array(self.matrix_to_global[:,2])
	
	def set_medium_index(self, new_index):
		self.q.set_medium_index(new_index)
		return self
	
	def get_medium_index(self):
		return self.q.medium_index
		
	def set_lambda(self, lambda0):
		self.q.lambda0=lambda0
		return self
	
	def get_lambda(self):
		return self.q.lambda0
		
	def free_drift(self, distance):
		self.q.drift(distance)
		self.total_drift+=distance
		self.x0+=self.direction()*distance
		self.total_drift_time+=distance*self.q.medium_index
		self.incoming_q=self.q.clone()
		self.footprint_q=self.q.clone()
		return self
		
	def localize(self, optic):			
		ar=Numeric.array
		dot=Numeric.dot
		tr=Numeric.transpose
		
		if hasattr(self,"local_x0"):
			raise "attempt to localize already localized beam at optic: "+optic.name
		self.incoming_direction=self.direction()
		cm=dot(optic.matrix_to_local, self.matrix_to_global)
		cm2=ar(cm[0:2,0:2])
		a,b,c,d=tuple(cm2.flat)
		cm2inv=ar(((d,-b),(-c,a)))/(a*d-b*c)
		self.incoming_q=self.q.clone()
		self.q.transform(cm2)
		self.footprint_q=self.q.clone() #keep a copy of the localized q for markers
		self.local_x0=dot(optic.matrix_to_local, self.x0-optic.center)
		self.local_direction=ar(cm[:,2])
		self.localize_tensor_transform=cm2inv
		self.local_polarization=dot(cm2, self.polarization)
		
	def globalize(self, optic):			
		ar=Numeric.array
		dot=Numeric.dot
		tr=Numeric.transpose
		
		cm=dot(optic.matrix_to_local, self.matrix_to_global)
		cm2=ar(cm[0:2,0:2])
		a,b,c,d=tuple(cm2.flat)
		cm2inv=ar(((d,-b),(-c,a)))/(a*d-b*c)
		self.q.transform(cm2inv)	
		self.polarization=dot(cm2inv, self.local_polarization)
		
		del self.local_x0, self.local_direction, self.localize_tensor_transform, self.local_polarization
		
	def transform(self, matrix):
		ar=Numeric.array
		dot=Numeric.dot
		self.matrix_to_global=dot(matrix, self.matrix_to_global)
			
	def update_q(self, optic):
		self.q.qw(optic)
			
	def __str__(self):
		return "beam: x0=%s, dir=%s, q=%s" % (Numeric.array_str(self.x0, precision=3, suppress_small=1),
				Numeric.array_str(self.direction(), precision=3, suppress_small=1),
				str(self.q))

	def __repr__(self):
		return self.__str__()
	
	def clone(self):
		return copy.deepcopy(self)
	
	def clone_no_marks(self):
		marks=self.marks
		del self.marks
		q=self.clone()
		self.marks=marks
		return q
	
	def shift_lambda(self, delta):
		self.set_lambda(self.get_lambda()+delta)
		return self
	
	def x(self, optic):
		optic.transform(self)
		return self #to make daisy-chains easy
	
	def mark(self, label=None): #record an entry in the marks array as to what we look like here
		if not self.marks.has_key(label):
			self.marks[label]=[] #empty list for this optic
			
		self.marks[label].append(self.clone_no_marks())
		self.marks["order"].append((label, len(self.marks[label])-1))

	def transformation_matrix_to_table(self, up):
		ar=Numeric.array
		dot=Numeric.dot
		tr=Numeric.transpose
		
		dir=self.matrix_to_global
		z=dir[:,2]
		x=cross(up,z)
		xmag=sqrt(dot(x,x))
		if xmag < 0.1 : #z is not even close to perpendicular to up!
			return 0, Numeric.identity(3).astype(Numeric.Float) #punt, return false and identity
		else:
			x=x/xmag
			y=cross(z,x) #may be tilted if the beam isn't completely level
			direction=tr(ar((x,y,z)))
			return 1, dot(tr(direction), dir) #coordinates are rational, return true and transform

	def transform_q_to_table(self, qi, up=(0,1,0)):
		dot=Numeric.dot
		tr=Numeric.transpose
		ok, matrix=self.transformation_matrix_to_table(up)
		if not ok : #z is not even close to perpendicular to up!
			theta, qix, qiy=qi.qi_moments()
			return qix, qiy #so just use principal moments	
		else:
			xyt=matrix[0:2,0:2]
			qiprime=dot(xyt, dot(qi.qit, tr(xyt)))
			return qiprime[0,0], qiprime[1,1]

class general_optic:
	
	def __init__(self, name, center=(0,0,0), angle=0, **extras):
		self.init(name, center, angle, extras)
		
	def init(self, name, center, angle, extras):
		
		self.abcd=Numeric.identity(2).astype(Numeric.Float) #predefine this as a default abcd matrix (in 2 dimensions)
		self.driftlength=0.0
		
		#copy any extra information directly into the class dictionary
		for i in extras.keys():
			self.__dict__[i]=extras[i]

		if len(center)==2: #2-d optics live in the x-z plane
			center=(center[0], 0, center[1])
		
		self.center=Numeric.array(center, Numeric.Float)
		self.name=name
		
		self.reset_angle(angle)
		self.post_init()
	
	def post_init(self):
		pass
	
	def add_info(self, **keys):
		for i in keys.keys():
			self.__dict__[i]=keys[i]
	
	def __str__(self):
		if self.name is not None:
			return (self.name+" "+Numeric.array_str(self.center, precision=3, suppress_small=1) +
					" "+Numeric.array_str(self.cosines, precision=3, suppress_small=1))
		else:
			return "Unnamed optic"
	
	def entrance_center(self):
		return self.center
	def exit_center(self):
		return self.center

	def transform_into_local(self, direction):
		return Numeric.dot(self.matrix_to_local, direction)
	
	def transform_into_global(self, local_direction):
		return Numeric.dot(self.matrix_to_global, local_direction)
		
	def check_hit_optic(self):
		"check to make sure local coordinates are within bounds of optic: default always true"
		return 1
	
	def localize_tensor(self, tens):
		dot=Numeric.dot
		tr=Numeric.transpose
		return dot(self.beam.localize_tensor_transform, dot(tens, tr(self.beam.localize_tensor_transform)))
	
	def globalize_transform(self, transform):
		dot=Numeric.dot
		return dot(self.matrix_to_global, dot(transform, self.matrix_to_local))
		
	def transform(self, beam, backwards=0):
		"transform(pos, direction, lambda) set up local coordinates, calls local_transform, and returns transform in global coordinates"
		self.beam=beam
		self.backwards=backwards #in case anyone needs to know (e.g. dielectric interfaces)
		beam.localize(self)
		self.check_hit_optic()
		self.local_transform()
		beam.globalize(self)
		del self.beam, self.backwards #get it out of our dictionary so it is an error if not localized
		return self
		
	def local_transform(self):
		#default is to use abcd_transform
		ar=Numeric.array
		dot=Numeric.dot
		tr=Numeric.transpose
		xp, yp, zp=tuple(self.beam.local_direction)  #it better be a small angle, so sin(theta)=theta for this to work
		dx, dy, dz=tuple(self.beam.local_x0) #dz should always be zero in these coordinates
		x1,xp1, y1, yp1 =self.abcd_transform((dx, xp, dy, yp))
		sx=xp1-xp #this is sort of a sine of a rotation angle of x about y
		cx=math.sqrt(1-sx*sx)
		sy=yp1-yp #this is sort of a sine of a rotation angle of y about x
		cy=math.sqrt(1-sy*sy)
		rot=dot(ar(((cx, 0, -sx),(0,1,0),(sx,0,cx))),ar(((1,0,0),(0,cy,-sy),(0,sy,cy))))
		self.beam.transform(self.globalize_transform(rot))
		self.q_transform()
		
	def intersect(self, from_point, from_direction):
		a,x  =  planesolve(from_point, from_direction, self.center, self.cosines)
		
		if a < 0:
			raise OpticDirectionError, "Optic found on backwards side of beam: "+str(self)+" incoming (x,y,z)="+\
			Numeric.array_str(from_point, precision=3, suppress_small=1)+" cosines="+\
			Numeric.array_str(from_direction, precision=3, suppress_small=1) + \
			(" distance = %.3f, point=(%.3f %.3f %.3f)" % (a, x[0], x[1], x[2]))
		return a, x

	def transport_to_here(self, beam):
		distance, x=self.intersect(beam.x0, beam.direction())
		beam.free_drift(distance)
		return self
		
	def format_geometry(self):
		return ("(x,y,z)=(%.3f, %.3f, %.3f) "%tuple(self.center))+str(self.euler)
	
	def reset_angle(self, angle):
		if type(angle) != types.TupleType: #convert conventional angle to rotation about y only
			self.euler=euler(angle, 0, 0)
			theta=angle
		else: #this should be a triple angles theta (yaw), eta (roll around absolute z), phi(roll around local z)			
			self.euler=euler(angle[0], angle[1], angle[2])
			theta=angle[0]

		self.matrix_to_global=self.euler.mat
		self.matrix_to_local=Numeric.transpose(self.euler.mat)
		self.cosines=self.transform_into_global(Numeric.array((0,0,1)))
	
	def old_rotate_in_local_frame(self, matrix_to_global):
		self.matrix_to_global=Numeric.dot(self.matrix_to_global, matrix_to_global)
		self.matrix_to_local=Numeric.transpose(self.matrix_to_global)
		self.cosines=self.transform_into_global(Numeric.array((0,0,1)))
		
	def rotate_in_local_frame(self, matrix_to_global):
		self.update_coordinates(self.center, self.center, self.globalize_transform(matrix_to_global))

	def update_coordinates(self, parent_center=None, reference_coordinate=None, matrix_to_global=None):
		if parent_center is None:
			parent_center=self.center
		if reference_coordinate is None:
			reference_coordinate=self.center
		if matrix_to_global is not None:				
			self.matrix_to_global=Numeric.dot(matrix_to_global, self.matrix_to_global)
			self.matrix_to_local=Numeric.transpose(self.matrix_to_global)
		self.center=parent_center+Numeric.dot(matrix_to_global, self.center-reference_coordinate)
		self.cosines=self.transform_into_global(Numeric.array((0,0,1)))
		self.euler=None
		return self
		
	def polygon(self):
		try:
			
			face_centered=hasattr(self,"face_centered")
				
			if hasattr(self,"justify"):
				justify=self.justify
			else:
				justify="center"
						
			if justify=="left":
				left=0
				right=-self.width
			elif justify=="right":
				left=self.width
				right=0
			else:
				try:
					# a fraction of 1 top justifies, 0 center justifies, and -1 bottom justifies
					fraction=(float(justify)+1)*0.5
					left=self.width*(1-fraction)
					right=-self.width*fraction
				except:
					left=self.width/2.0
					right=-left
			
			if hasattr(self,"height"):
				height=self.height
			else:
				height=self.width	
							
			top=height/2.0
			bottom=-height/2.0
			
			if(face_centered):
				front=0
				back=self.thickness
			else:
				front=-self.thickness/2.0
				back=-front
							
			baserect=Numeric.array((
					(left, right, right, left, left, left, right, right, left, left, right, right), 
					(top, top, bottom, bottom, top, top, top, bottom, bottom, top, bottom, bottom),
					(front, front, front, front, front, back, back, back, back, back, back, front)))
			return Numeric.transpose(Numeric.dot(self.matrix_to_global, baserect))+self.center
			
		except:
			return None
	
	def polygon_list(self):
		return [self.polygon()] #for simple optics, the list has one item
	
	def place_between(self, from_obj, to_obj, distance):
		"place_between(from, to, distance, set_info) puts the optic between object or coordinate 'from' and 'to'"
		if isinstance(from_obj, general_optic):
			fc=from_obj.exit_center()
			fn=from_obj.name
		else:
			fc=Numeric.array(from_obj) #assume is is a tuple or array
			fn=Numeric.array_str(fc, precision=3, suppress_small=1)
		if isinstance(to_obj, general_optic):
			tc=to_obj.entrance_center()
			tn=to_obj.name
		else:
			tc=Numeric.array(to_obj) #assume is is a tuple or array
			tn=Numeric.array_str(tc, precision=3, suppress_small=1)
		
		dx1=tc-fc
		dx1hat=dx1/vec_mag(dx1)
		if distance > 0.0:
			self.center=fc+distance*dx1hat
		else:
			self.center=tc+distance*dx1hat
			
		x,y,z=tuple(dx1.flat)
		rho=math.sqrt(x*x+y*y)
		theta=math.atan2(rho,z)/deg
		eta=math.atan2(y,x)/deg
		self.reset_angle((theta, eta, 0))		
		return self #make daisy-chaining easy!

	def set_direction(self, from_obj, to_obj):
		"set_direction(from, to, set_info) points the mirror from object or coordinate 'from' to 'to'"
		if isinstance(from_obj, general_optic):
			fc=from_obj.exit_center()
			fn=from_obj.name
		else:
			fc=Numeric.array(from_obj) #assume is is a tuple or array
			fn=Numeric.array_str(fc, precision=3, suppress_small=1)
		if isinstance(to_obj, general_optic):
			tc=to_obj.entrance_center()
			tn=to_obj.name
		else:
			tc=Numeric.array(to_obj) #assume is is a tuple or array
			tn=Numeric.array_str(tc, precision=3, suppress_small=1)
		
		self.looking_from_name=fn
		self.looking_to_name=tn
		self.looking_from_obj=from_obj
		self.looking_to_obj=to_obj
		
		dx1=fc-self.center
		dx2=tc-self.center
		dx1hat=dx1/vec_mag(dx1)
		dx2hat=dx2/vec_mag(dx2)
		halfway=(dx1hat+dx2hat)/vec_mag(dx1hat+dx2hat)
		x,y,z=tuple(halfway.flat)
		rho=math.sqrt(x*x+y*y)
		theta=math.atan2(rho,z)/deg
		eta=math.atan2(y,x)/deg
		self.reset_angle((180+theta, eta, 0))	
		self.perp=Numeric.array((dx1hat[1]*dx2hat[2]-dx1hat[2]*dx2hat[1],
				dx1hat[2]*dx2hat[0]-dx1hat[0]*dx2hat[2],
				dx1hat[0]*dx2hat[1]-dx1hat[1]*dx2hat[0]))  #cross product perpendicular to plane of incidence, for reference
		return self #make daisy-chaining easy

	def rotate_axis(self, axis_theta):
		self.rotate_in_local_frame(matrix_to_global=eulermat(0.0, 0.0, axis_theta))
		return self

	def tilt_off_axis(self, angle):
		if type(angle) is not types.TupleType:
			theta=angle
			eta=0.0
		else:
			theta=angle[0]
			eta=angle[1]
		self.rotate_in_local_frame(matrix_to_global=eulermat(theta, eta, 0.0))
		return self
		
	def clone(self, newname=None):
		a=copy.deepcopy(self)
		if newname is not None:
			a.name=newname
		return a
		
	def format_name(self):
		if self.name is None:
			return ""
		else:
			return self.name + ": "
	
class base_reflector:
	"reflector is a class for mirror-like objects... the default is a mirror"
	def local_transform(self):
		self.beam.transform(self.globalize_transform(Numeric.array(((1.0,0,0),(0,1.0,0),(0,0,-1.0)))))
	
	def post_init(self):
		self.reflector_info=self.format_geometry()

class base_lens:
	def post_init(self):
		f=self.f
		if type(f)==types.TupleType:
			f1, f2 = f
			if f1 is not None:
				d1=-1.0/f1
			else:
				d1=0.0
			if f2 is not None:
				d2=-1.0/f2
			else:
				d2=0.0
		else:
			if f is not Infinity:
				d1=d2=-1.0/f
			else:
				d1=d2=0.0
		
		self.d1=d1
		self.d2=d2
		self.strength=Numeric.array(((d1,0),(0,d2)))
		self.info_str()
		
	def mat2(self):
		if self.d1==self.d2:
			return Numeric.array(((1.0,0.),(-1.0/self.d1,1)))
		else:
			raise "Attempt to get 2x2 matrix of non-spherical lens" 

	def mat4(self):
		mxx, mxy, myx, myy=[Numeric.array(((1.,0.),(d,1.))) for d in self.strength.flat]
		return Numeric.array(((mxx,mxy),(myx,myy)))
	
	def rotate_axis(self, axis_theta):
		general_optic.rotate_axis(self, axis_theta)
		self.info_str()
		return self
		
	def q_transform(self):
		self.beam.q.focus(self.strength) #it's easy!
	
	def abcd_transform(self, vec):
		x, xp, y, yp = tuple(vec)
		
		dxx,dxy,dyx,dyy=tuple(self.localize_tensor(self.strength).flat)
		return x, xp+dxx*x+dxy*y, y, yp+dyy*y+dyx*x
	
	def info_str(self): #short format string
		if self.d1==self.d2:
			if self.d1 != 0:
				self.base_lens_info=("f=%.1f"%(-1.0/self.d1))
			else:
				self.base_lens_info=("f= Infinity")

		else:
			if self.d1==0:
				f1=""
			else:
				f1="f1=%.1f " % (-1.0/self.d1)
			if self.d2==0:
				f2=""
			else:
				f2="f2=%.1f " % (-1.0/self.d2)
			self.base_lens_info=f1+f2+("axis theta=%.1f"%(self.axis_theta))

	def __str__(self):
		return self.format_name()+self.base_lens_info

class reflector(general_optic, base_reflector):
	"reflector is a class for mirror-like objects... the default is a mirror"
		
	def post_init(self):
		if not hasattr(self,"face_centered"):
			self.face_centered=1
		if not hasattr(self,"width"):
			self.width=0.0254
		if not hasattr(self,"height"):
			self.height=self.width
		if not hasattr(self,"thickness"):
			self.thickness=0.01
		base_reflector.post_init(self)
		
	def __str__(self):
		return self.format_name()+self.reflector_info
		
	def local_transform(self):
		base_reflector.local_transform(self)
	
	def q_transform(self):
		pass
		
	def set_direction(self, from_obj, to_obj, set_info=1):
		"set_direction(from, to, set_info) points the mirror from object or coordinate 'from' to 'to'"
		general_optic.set_direction(self, from_obj, to_obj)
		
		if set_info:
			self.reflector_info=("looking from %s to %s " % 
					(self.looking_from_name, self.looking_to_name))+str(self.euler)
		return self #make daisy-chaining easy
		
class null_optic(general_optic):
	def post_init(self):
		if not hasattr(self,"width"):
			self.width=0.0254
		if not hasattr(self,"thickness"):
			self.thickness=0.01
			
	def abcd_transform(self, vec):
		return vec

	def q_transform(self):
		pass
		
class lens(general_optic, base_lens):
	def post_init(self):
		base_lens.post_init(self)
		#default lens draws as 1" diameter x 5 mm thick
		if not hasattr(self,'thickness'): self.thickness=0.005
		if not hasattr(self,'width'): self.width=0.0254
		if not hasattr(self,'height'): self.height=self.width

	def __str__(self):
		return base_lens.__str__(self)+" "+self.format_geometry()	

class dielectric_interface(general_optic, base_lens):
	""" warning... as of 9/24/2002, this is not even close to done. Do not use"""
	
	def post_init(self):
		#raise "Warning!", "dielectric_interface is not ready to use! Do not try this."
		self.last_lambda=0.0 
		if not hasattr(self, "external_index"):
			self.external_index=get_ambient_index
		
	def local_transform(self):
		
		dot=Numeric.dot
		
		self.update_strength()
		
		kx, ky, kz=self.beam.local_direction
		
		#must rotate by the appropriate refraction angle _in the plane of incidence_ 
		sin_theta_in=math.sqrt(kx*kx+ky*ky) #incoming k had better be a unit vector!
		if abs(sin_theta_in) < 1e-6: return
		cphi, sphi = kx / sin_theta_in, ky/sin_theta_in
		
		n1=self.ambient_index
		n2=self.medium_index
				
		if kz < 0.0:
			n1, n2 = n2, n1 #we're going backwards through the interface
			sin_theta_in=-sin_theta_in
		
		self.final_index=n2
		
		sin_theta_out= n1/n2*sin_theta_in
		dtheta=math.asin(sin_theta_out)-math.asin(sin_theta_in)
		c,s = math.cos(dtheta), math.sin(dtheta)
		phimat= Numeric.array(((cphi,-sphi, 0),( sphi, cphi, 0),(0,0,1.)))
		phimati=Numeric.array(((cphi, sphi, 0),(-sphi, cphi, 0),(0,0,1.)))
		
		#this matrix is a rotation by by phi to put the incident beam on x, rotate about x-z by dtheta, then undo the phi
		self.beam.transform(self.globalize_transform(dot(phimat,dot(Numeric.array(((c,0,s),(0,1,0),(-s,0,c))), phimati))))
		self.q_transform()
		

	def update_strength(self):
		oldlambda=self.beam.get_lambda()
		if oldlambda==self.last_lambda: return
		
		if type(self.external_index) is types.FunctionType:
			self.ambient_index=self.external_index(oldlambda)
		else:
			self.ambient_index=self.external_index
		n1=self.ambient_index
		
		if type(self.index) is types.FunctionType:
			self.medium_index=self.index(oldlambda)
		else:
			self.medium_index=self.index
			
		n2=self.medium_index
		
		if hasattr(self, "radius"):
			raise "Warning...", "Curved dielectric interfaces don't work yet! try again later"
			if type(self.radius) is not types.TupleType:
				rx=ry=self.radius
			else:
				rx,ry=self.radius
			if rx is not Infinity:
				cx=(n1-n2)/(rx*n2)
			else:
				cx=0.0
			if ry is not Infinity:
				cy=(n1-n2)/(ry*n2)
			else:
				cy=0.0
			self.d_scalar=n1/n2	#always diagonal, for scalar media, so we don't have to localize 	
			self.strength=Numeric.array(((cx,0),(0,cy)))
		else:
			self.strength=Numeric.array(((0.,0.),(0.,0.)))
			
		self.last_lambda=oldlambda
		
	def q_transform(self):
		self.beam.set_medium_index(self.final_index)
		kx, ky, kz=self.beam.local_direction
		if kz > 0.0:
			self.beam.q.focus(self.strength)
		else: #going the wrong way!
			self.beam.q.focus(self.strength*(-self.medium_index/self.ambient_index))
			

class spherical_mirror(reflector, base_lens):
	def post_init(self):
		base_reflector.post_init(self)
		base_lens.post_init(self)
		#default lens draws as 1" diameter x 5 mm thick
		if not hasattr(self,'thickness'): self.thickness=0.005
		if not hasattr(self,'width'): self.width=0.0254
		if not hasattr(self,'height'): self.height=0.0254
		
	def local_transform(self, beam):
		general_optic.local_transform(self) #pick up our abcd transform for the lens
		reflector.local_transform(self) #and do the reflection
	
	def __str__(self):
		return self.format_name()+self.reflector_info+" "+self.base_lens_info+" "+self.format_geometry()
if 0:
	print "\n\n****start mirror tests"
	mir1=reflector("reflector1", center=(0,0,1))
	mir2=reflector("reflector2", center=(0,1,1))
	mir1.set_direction((0,0,0), mir2)
	mir2.set_direction(mir1, (1,1,1))
	mybeam=beam((0,0,0), qtens(1.054e-6, r=Infinity, w=.002))
	mybeam.q.qit[0,0]*=0.25 
	for optic in (mir1, mir2):
		print optic.name
		optic.transport_to_here(mybeam)
		print mybeam
		optic.transform(mybeam)
		print mybeam

if 0:
	print "\n\n****start lens tests"
	optic=lens("reflector", center=(0,0,1), f=0.25)
	mybeam=beam((0,0,0), qtens(1.054e-6, r=Infinity, w=.002)) 
	print mybeam
	optic.transport_to_here(mybeam)
	print mybeam
	optic.transform(mybeam)
	print mybeam


class grating(reflector):
	"grating is a reflector which diffracts. The ruling is along the y axis"
	
	def local_transform(self):
		
		kx, ky, kz=self.beam.local_direction
		
		theta=math.atan2(kx, kz)
		littrow, out=self.angles(theta, self.beam.get_lambda())
		
		#remember, outgoing angle is _not_ outgoing direction... beam has changed from incoming to outgoing, too
		dtheta=out+theta
		s,c=math.sin(dtheta), math.cos(dtheta)
		
		#this funny matrix is a rotation by dtheta, followed by a reflection in z
		self.beam.transform(self.globalize_transform(Numeric.array(((c,0,-s),(0,1,0),(-s,0,-c)))))
						
	def degree_angles(self, theta, lam):
		#return angle info for incoming and outgoing angles in degrees (for user convenience)
		litt, beta = self.angles(theta*deg, lam)
		return litt/deg, beta/deg
		
	def angles(self, theta, lam):
		#return angle info for incoming and outgoing angles in degrees (for user convenience)
		if theta<0:
			sgn=-1
		else:
			sgn=1
		gratparm=self.order*lam*self.pitch*sgn
		littrow=asin(gratparm/2.0)
		phi=(((theta + math.pi/2.0) % math.pi ) - math.pi/2.0) #always wrap to front of grating!
		beta0=asin(gratparm-sin(phi))
		self.parafocus_scale=math.cos(phi)/math.cos(beta0) #angular magnification along x
		#print "grating angles: ", lam, theta/deg, littrow/deg, beta0/deg
		return (littrow, beta0)
		
	def __str__(self):
		return reflector.__str__(self)+("pitch= %.3f" % self.pitch)

class key_tag:
	def __init__(self, key):
		self.key=key

class backwards(key_tag):
	pass

def get_tagged_key(object):
	if isinstance(object, key_tag):
		return object.key
	else:
		return object

class composite_optic(general_optic):
	"composite_optic is a container for a list of general_optics which can be rotated, etc. as a group"
	
	def __init__(self, name, optics_dict, optics_order, reference_coordinate=None,
			center=(0,0,0), angle=0.0, **extras):
		self.init( name, optics_dict, optics_order, reference_coordinate, center, angle, extras)
		
	def init(self, name, optics_dict, optics_order, reference_coordinate, center, angle, extras):
		
		if reference_coordinate is None: #refer all coordinates to first optic in list by default
			reference_coordinate=optics_dict[optics_order[0]].entrance_center()
			
		self.optics_dict={}
		self.optics_order=list(optics_order)

		for k in optics_order: #only pull out keys which are used, and let python eliminate duplicates
			self.optics_dict[get_tagged_key(k)]=optics_dict[get_tagged_key(k)]
						
		general_optic.init(self, name, center=center, angle=angle, extras=extras)
		for k in self.optics_dict.keys():
			self.optics_dict[k].update_coordinates(center, reference_coordinate, self.matrix_to_global)
	
	def exit_center(self):
		return self.optics_dict[self.exit_optics_tags()[1]].center
	
	def entrance_center(self):
		return self.optics_dict[self.entrance_optics_tags()[0]].center
				
	def update_coordinates(self, new_center, parent_reference, matrix_to_global):
		for k in self.optics_dict.keys():
			self.optics_dict[k].update_coordinates(new_center, parent_reference, matrix_to_global)
		general_optic.update_coordinates(self, new_center, parent_reference, matrix_to_global) #update our top-level matrices, etc.
		
	def mark_label(self, opticname):
		return (self,opticname)
				
	def transform(self, beam, back=0):
		if back:
			order=copy.copy(self.optics_order)
			order.reverse()
		else:
			order=self.optics_order
		
		for opticname in order:
			optic=self.optics_dict[get_tagged_key(opticname)]
			if not isinstance(optic, composite_optic): #only mark primitive optics
				waist=beam.q.next_waist()
				current_drift=beam.total_drift
				optic.transport_to_here(beam)
				distance=beam.total_drift-current_drift
				if waist > 0 and waist < distance:
					beam.free_drift(waist-distance)
					beam.mark("waist")
					optic.transport_to_here(beam)
			if isinstance(opticname, backwards):
				optic.transform(beam, not back) #'backwards' must toggle the reversed flag, so backwards-backwards is forwards
			else:
				optic.transform(beam, back)
				
			if not isinstance(optic, composite_optic): #only mark primitive optics
				beam.mark(self.mark_label(get_tagged_key(opticname)))
			
	def polygon_list(self):
		l=[]
		for o in self.optics_dict.keys():
			#use additions semantics to flatten nesting of lists
			l=l+self.optics_dict[o].polygon_list()
		return [poly for poly in l if poly is not None]	
	
	def exit_optics_tags(self):
		return self.optics_order[-2], self.optics_order[-1]
	
	def entrance_optics_tags(self):
		return self.optics_order[0], self.optics_order[1]

	def set_exit_direction(self, look_to):
		l=self.optics_dict
		e1, e2 = self.exit_optics_tags()
		l[e2].set_direction(l[e1], look_to)

	def set_entrance_direction(self, look_to):
		l=self.optics_dict
		e1, e2 = self.entrance_optics_tags()
		l[e1].set_direction(look_to, l[e2])

	def rotate_to_axis(self, from_obj):
		if isinstance(from_obj, general_optic):
			fc=from_obj.exit_center()
			fn=from_obj.name
		else:
			fc=Numeric.array(from_obj) #assume is is a tuple or array
			fn=Numeric.array_str(fc, precision=3, suppress_small=1)
		
		dx1=self.center-fc
		x,y,z=tuple(dx1)
		rho=math.sqrt(x*x+y*y)
		theta=math.atan2(rho,z)/deg
		eta=math.atan2(y,x)/deg
		if abs(eta)>90: #prefer signed theta and unsigned eta
			eta=180-eta
			theta=-theta
		self.update_coordinates(self.center, self.center, euler(theta, eta, 0).mat)	
		return self #make daisy-chaining easy

	def __getitem__(self, tag):
		if type(tag) is types.TupleType and tag[0] is self:
			return self.optics_dict[tag[1]] #in case the tag had our object identifier attached
		else:
			return self.optics_dict[tag]

	def __len__(self):
		return len(self.optics_order)
		
if 0:
	print "testing grating"
	optic=grating("grating", pitch=1.5e6, angle=38.0, order=1)
	mybeam=beam((0,0,0), qtens(1.054e-6, r=Infinity, w=.002)) 
	print mybeam
	optic.transport_to_here(mybeam)
	print mybeam
	optic.transform(mybeam)
	print mybeam
	
class optics_trace:
	def __init__(self, marks, color=None, **extras):
		self.color=color
		self.marks=marks
		#copy any extra information directly into the class dictionary
		for i in extras.keys():
			self.__dict__[i]=extras[i]
	
	def __getitem__(self, index):
		if type(index) is types.IntType:
			tag, n = self.marks["order"][index]
			return self.marks[tag][n]
		else:
			return self.marks[index[0]][index[1]]

	def __len__(self):
		return len(self.marks["order"])
	
def trace_path(optics_system, beam):
	beam.total_drift=0.0		
	try:
		optics_system.transform(beam)	
	except:
		traceback.print_exc()
		pass #just return an incomplete trace if something goes wrong!
	return optics_trace(beam.marks)

class phase_plate(null_optic):
	def post_init(self):
		theta=self.phase_angle*math.pi/180.0
		self.x_phase_shift=complex(math.cos(theta), math.sin(theta))
	
	def local_transform(self):
		self.beam.local_polarization[1]*=self.x_phase_shift
	
	def info_str(self):
		return "Phase plate: x shift = %.1f" % self.phase_angle
		
class halfwave_plate(phase_plate):
	def post_init(self):
		self.phase_angle=180.0
		phase_plate.post_init(self)

class quarterwave_plate(phase_plate):
	def post_init(self):
		self.phase_angle=90
		phase_plate.post_init(self)
	
class faraday_rotator(null_optic):
	def post_init(self):
		pass

	def info_str(self):
		return "Faraday Rotator: rotation = %.1f" % self.rotation

	def local_transform(self):
		theta=self.rotation*self.beam.local_direction[2]
		s,c=sincosdeg(theta)
		#if beam is going forward, shift (+), if backwards, shift (-)
		#where forward/backwards comes from the local z-component of the beam direction
 		self.beam.local_polarization=Numeric.dot(Numeric.array(((c,-s),(s,c))), self.beam.local_polarization)

class dielectric_trapezoid(composite_optic):
	ENTRANCE='entrance'
	CENTER='center'
	EXIT='exit'
	def __init__(self, name, entrance_face_angle=0.0, exit_face_angle=0.0, thickness=0.05, width=0.025, index=None, center=(0.,0.,0.), angle=0.0, **extras):
		
		if index is None:
			raise "No index of refraction (scalar or function) given as index=foo creating dielectric trapezoid"
		
		my=dielectric_trapezoid
		self.thickness=thickness
		self.width=width
		
		self.front_face_offset=math.tan(entrance_face_angle*deg)*width*0.5
		self.back_face_offset=math.tan(exit_face_angle*deg)*width*0.5
		
		optics={}
		optics[my.ENTRANCE]=dielectric_interface("entrance face", angle=entrance_face_angle, index=index,
				center=(0,0,-0.5*thickness))
		optics[my.CENTER]=null_optic("trapezoid center", center=(0,0,0), width=width*0.05, thickness=width*0.05)
		optics[my.EXIT]=dielectric_interface("exit face", angle=exit_face_angle+180.0, index=index,
				center=(0,0,0.5*thickness) )
		order=[my.ENTRANCE, my.CENTER, my.EXIT]
		
		#specify reference as (0,0,0) to make our coordinates referenced to (0,0,0)
		composite_optic.init(self, name, optics, order, (0,0,0), center, angle, extras=extras )
		
	def polygon(self):
		if hasattr(self,"height"):
			height=self.height
		else:
			height=self.width	
						
		top=height/2.0
		bottom=-height/2.0
		
		front=-self.thickness/2.0
		back=-front
		
		fo=self.front_face_offset
		bo=self.back_face_offset
		
		left=0.5*self.width
		right=-left
		
		baserect=Numeric.array((
				(left, right, right, left, left, left, right, right, left, left, right, right), 
				(top, top, bottom, bottom, top, top, top, bottom, bottom, top, bottom, bottom),
				(front-fo, front+fo, front+fo, front-fo, front-fo, back-bo, back+bo, back+bo, back-bo, back-bo, back+bo, front+fo)))
		return Numeric.transpose(Numeric.dot(self.matrix_to_global, baserect))+self.center

	def polygon_list(self):
		my=dielectric_trapezoid
		return [self.polygon(), self.optics_dict[my.CENTER].polygon()] 


if 0:
	print "\n\n****start phase plate tests"
	optic=halfwave_plate("phase", center=(0,0,1)).rotate_axis(30)
	mybeam=beam((0,0,0), qtens(1.054e-6, r=Infinity, w=.002)) 
	print mybeam
	optic.transport_to_here(mybeam)
	print mybeam, mybeam.polarization
	optic.transform(mybeam)
	print mybeam, mybeam.polarization
