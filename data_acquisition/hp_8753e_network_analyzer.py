"A basic wrapper for some of the more useful external fuctionality for the Hewlett-Packard / Agilent 8753E-class Network Analyzers"
_rcsid="$Id: hp_8753e_network_analyzer.py,v 1.1 2003-08-06 18:22:29 mendenhall Exp $"

from vxi_11 import vxi_11_connection
import struct
import Numeric
import time
import array

_bigendian=(struct.unpack('H','\00\01')[0] == 1)

class HP8753_Error(Exception):
	pass

class hp_8753e(vxi_11_connection):
	
	idn_head="HEWLETT PACKARD,8753E" 	
				
	default_lock_timeout=10000

	def check_errors(self):
		errorlist=[]
		while(1):
			err, stat=self.read_status_byte()
			has_error=stat & 8
			if not has_error: break
			err, reason, result=self.transaction("OUTPERRO;\n")
			s1, s2 = result.strip().split(',', 1)
			errorlist+=[(int(s1),s2[1:-1])]
		if errorlist:
			raise HP8753_Error(errorlist)
		
	def setup_simple_sweep(self, min=1e4, max=3e9, count=1601, sweeptype='lin', power=10.0):
		if sweeptype=='lin': #when switching frequencies, always disable transform mode and corrections
			s="TIMDTRAN OFF;CORR OFF;LINFREQ;"
		else:
			s="TIMDTRAN OFF;CORR OFF;LOGFREQ;"
		s+="POIN %d;STAR %.8e Hz;STOP %.8e Hz;POWE %.2f"%(count, min, max, power)
		self.lock()
		try:
			self.write(s)
			self.check_errors()
		finally:
			self.unlock()
	
	def enable_transform(self, tmin=-1e-8, tmax=1e-8, mode='impulse'):
		if mode=='impulse':
			s='LOWPIMPU'
		elif mode=='step':
			s='LOWPSTEP'
		elif mode=='bandpass':
			s='BANDPASS'
		else:
			raise HP8753_Error("Bad transform mode: "+mode)
		
		self.write("SETF;TIMDTRAN ON;%s;REAL;STAR %.8e;STOP %.8e;"%(s, tmin, tmax))
		self.check_errors() #make sure it worked
	
	def disable_transform(self):
		self.write("TIMDTRAN OFF;")

	def recall_setup(self, setup_reg):
		self.write("RECAREG %d;" % setup_reg)
		self.check_errors()
	
	def get_converted_data(self, command):
		"return array of  data from active channel as a native-endian double array"
		err, reason, result=self.transaction("FORM3;"+command)
		if result[:2] != "#A":
			raise HP8753_Error("Bad data string returned: "+repr(result))
		count=struct.unpack(">H",result[2:4])[0]
		if len(result) != count+4:
			raise HP8753_Error("Bad data string length:  expected %d, got %d"%(count+4, len(result)))
		values=array.array('d')
		values.fromstring(result[4:])
		if not _bigendian:
			values.byteswap()		
		return values
		
	def get_complex_data(self):
		"return array of complex data from active channel"
		values=self.get_converted_data("OUTPDATF;")
		cv=Numeric.fromstring(values.tostring(), Numeric.Complex64)
		return cv

	def get_real_display_data(self):
		"return array of real data from active channel (assuming channel is a real format)"
		values=self.get_converted_data("OUTPFORF;")
		cv=Numeric.fromstring(values.tostring(), Numeric.Float64)
		return cv

	def get_complex_display_data(self):
		"return array of complex data from active channel"
		values=self.get_converted_data("OUTPFORF;")
		cv=Numeric.fromstring(values.tostring(), Numeric.Complex64)
		return cv

	def generate_xaxis_data(self):
		"figure out the correct x axis for the current data set"
		self.lock()
		try:
			info=[self.transaction(command)[2] for command in ["POIN?;","STAR?;","STOP?;","LINFREQ?;","LOGFREQ?;","TIMDTRAN?;"] ]
		finally:
			self.unlock()
		linf, logf, transform = [int(info[i]) for i in [3,4,5]] 
		fpoints, start, stop = [float(info[i]) for i in [0, 1,2]]
		points=int(fpoints)
		
		basearray=Numeric.array(range(points), Numeric.Float)*(1.0/(points-1)) #array on [0,1]
		
		if transform or linf:
			return basearray*(stop-start)+start
		else:
			return Numeric.exp(basearray*Numeric.log(stop/start))*start
							
if __name__=='__main__':
	from pyx import *
	import os
	import sys
	
	class graphxy(graph.graphxy):
		" a local pyx.graph.graphxy with some extra housekeeping and a spare copy of LaTeX"
		def __init__(self, **kwargs):
			graph.graphxy.__init__(self,**kwargs)
			self.settexrunner(text.texrunner(mode='tex'))
			self.latex=text.texrunner(mode='latex')
			
		def display(self, file):
			self.writetofile(file)
			os.system('rm %s.* %s.*' % (self.texrunner.texfilename, self.latex.texfilename) )
					
			if sys.platform=='darwin': #use pdf display with open on Mac
				os.system("epstopdf %s.eps;rm %s.eps;open %s.pdf"% (3*(file,)))
			else:
				os.system("gv %s.eps &"%file)
	
		def latex_text(self, *args):
			self.insert(self.latex.text(*args))
		
			
	analyzer=hp_8753e(host="***REMOVED***", 
		 device="gpib0,3",  timeout=4000, device_name="network_analyzer", raise_on_err=1)
	print analyzer.idn	
	analyzer.lock()
	try:
		analyzer.setup_simple_sweep(min=1e6,max=6e9, sweeptype='lin', power=10.0)
		analyzer.write("S11;")
		analyzer.enable_transform(tmin=-5e-9, tmax=15e-9, mode='impulse')
		analyzer.check_errors()
		cv=analyzer.get_real_display_data()
		
		xaxis=analyzer.generate_xaxis_data()*1e9
		g=graphxy(width=30,x=graph.linaxis(title='time (ns)', min=-5, max=15), y=graph.linaxis())
		g.plot(graph.data(data.data(zip(xaxis, cv)),x=0,y=1), graph.line())
		g.texrunner.settex(lfs='foils17pt')
		g.display('test')
		
		
		
	finally:
		analyzer.unlock()
		del analyzer
	
