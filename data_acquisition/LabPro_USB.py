"LabPro_USB supports connections of the Vernier LabPro system via USB"

_rcsid="$Id: LabPro_USB.py,v 1.1 2003-06-04 17:26:50 mendenhall Exp $

import LabPro
LabPro import LabPro, LabProError, _bigendian

import array
import Numeric

import time
import sys

import os
import threading
import traceback
import fcntl
import struct

class USB_data_mixin:
	"mixin class for RawLabPro to process USB-blocked data from LabPro"
				
	def get_data_binary(self, chan=0, points=None, scaled_range=None):
		"get specified channel data, either as raw integers if scaled_range=None, or as floats scaled to specified range in Numeric array"
		#note... in USB mode, LabPro doesn't use checksums, so this is modified from the usual
		if points is None:
			points=self.get_system_config()['dataend']
		self.command('data_control',chan,0,0,0,1)
		self.send_string('g')
		chunklen=(2*points+63) & (-64) #round up to next 64-byte block, since LabPro always returns 64-byte multiples in USB
		s=''
		empties=0
		while(empties<5 and len(s)<chunklen):
			time.sleep(0.02)
			newdata=self.read(chunklen-len(s), mode=1) #mode=1 stops stripping of nulls from end of 64 bytes packets
			s+=newdata
			if newdata:
				empties=0
			else:
				empties+=1
		
		if empties:
			raise LabProTimeout('timeout getting binary data, got %d bytes, expected %d, stuff looks like: %s' %(len(s), chunklen, repr(s[:100])))		
			
		data=array.array('H')
		data.fromstring(s[:2*points])
		if not _bigendian:
			data.byteswap()
		if scaled_range is None:
			return data #return raw integers	
		else:
			return self.scale_binary_data(data, scaled_range) 

	def binary_mode(self, blocking_factor=4):
		"set how many readings the LabPro returns in a single USB realtime block, and transfer in binary"
		self.binary_blocking_factor=blocking_factor
		self.command(4,0,-1,blocking_factor) #blocked binary data
		self.__saved_realtime_fragment=''
		
	def parse_binary(self, data, channels):
		"convert LabPro realtime binary string, with USB data blocking, to integer list"
		#print repr(data)
		format='>'+channels*'H'+'L' #data+time
		chunklen=2*channels+4
		l=[]
		for block in range(len(data)//64):
			for chunk in range(self.binary_blocking_factor):
				offset=block*64+chunk*16
				l.append(struct.unpack(format,data[offset:offset+chunklen]))
		return l
		
	def get_data_binary_realtime(self, channels=1):
		"read the most recently collected realtime data from the LabPro and return as a list, using binary transfers and no scaling"
		s=self.__saved_realtime_fragment
		while(len(s)<64):
			time.sleep(0.005)
			s=s+self.read(maxlen=None, mode=1) #binary, don't chop zeros
			if not s:
				return [] #no data if string is still blank
		leftovers=len(s) % 64
		if leftovers:
			self.__saved_realtime_fragment=s[-leftovers:]
			s=s[:-leftovers]
		else:
			self.__saved_realtime_fragment=''
	
		return self.parse_binary(s, channels)

class USB_Mac_mixin:
	"mixin class for RawLabPro to allow operation of LabPro via USB port on Macintosh OSX using pipe server"
		
	def setup_serial(self,port_name=None):
		self.usb_send, self.usb_recv, self.usb_err=os.popen3(os.path.join(os.path.dirname(__file__),"LabProUSB"),'b',0)
		
		fcntl.fcntl(self.usb_recv, fcntl.F_SETFL, os.O_NONBLOCK) #pipes must be nonblocking
		fcntl.fcntl(self.usb_err, fcntl.F_SETFL, os.O_NONBLOCK) #pipes must be nonblocking
		self.__keep_running=1
		self.status_monitor=threading.Thread(target=self.read_status, name='USB LabPro status monitor')
		self.status_monitor.start()
		self.__saved_realtime_fragment=''
		
	def check_usb_status(self):
		if not self.status_monitor.isAlive():
			raise LabProError("LabPro USB server has died...")
			
	def high_speed_serial(self):
		"not implemented on USB"
		return 

	def high_speed_setup(self):
		"not implemented onUSB"
		return 

	def set_port_params(self):
		"not implemented on USB"
		return 

	def read(self, maxlen=None, mode=None):
		"read data from USB.  If mode is None or 0, strip trailing nulls for ASCII, otherwise leave alone"
		self.check_usb_status()
		try:
			if maxlen is None:
				res=self.usb_recv.read() 
			else:
				res=self.usb_recv.read(maxlen)
			if not mode:
				zp=res.find('\0')
				if zp>=0:
					res=res[:zp] #trim any nulls
			return res
			
		except IOError:
			if sys.exc_value[0]==35: #this error is returned on a nonblocking empty read
				return '' #just return empty data
			else:
				raise LabProError("USB server disconnected unexpectedly")

	def read_status(self):
		"monitor the pipe server's stderr in a thread"
		while(self.__keep_running):
			try:
				res=self.usb_err.read()
				if res: print >> sys.stderr, "LabPro USB message: ", res
				if res.find("****EXITED****") >= 0: break #if this thread dies, LabPro has gone away somehow
			except IOError:
				if sys.exc_value[0]==35: #this error is returned on a nonblocking empty read
					time.sleep(1) #so just wait and try again
				else:
					raise
		
	def write(self, data):
		self.check_usb_status()
		self.usb_send.write(data)
		time.sleep(0.1) #give a little extra time for message passing
		
	def close(self):
		self.__keep_running=0
		try:
			self.usb_send.write("****QUIT****\n")
		except:
			pass
	
class LabPro_Mac_USB(USB_data_mixin, USB_Mac_mixin, LabPro):
	def __init__(self):
		LabPro.__init__(self,'')

if __name__=='__main__':
	
	class LabProData(Exception):
		pass
	
	def monitor_humidity(filepath="/Users/marcus/Public/humidity.txt", sample_period=30.0):
		
		badloops=0
		lp=LabPro_Mac_USB()
		
		try:	
			lp.wake()
			lp.reset()
			print lp.get_system_config() #just to vacuum the connection
			lp.setup_channel(chan=2, operation=14, equation=0) #set up 0-5 volts readback for humidity sensor
			lp.command(107, .001, 1000) #setup 1000-1 oversampling for 1 second		
			while(1):
				time.sleep(sample_period)
				try:
	
					lp.setup_data_collection(samptime=2.0, numpoints=1, filter=0)
					lp.wait_for_data_done(flasher=0)
					lp.binary_mode()
					value=lp.get_data_binary(chan=2, scaled_range=(0.,5.))						
					if len(value)>1:
						raise LabProData("got too much data...: "+str(value))
					value=value[0]
		
					grabtime=time.time()
					try:
						if filepath is not None:
							f=file(filepath,"a")
						else:
							f=sys.stdout
							
						print >> f, "%s\t%.1f\t%.3f\t%.1f" % (time.asctime(), grabtime, value, value*31.47-31.74)
						
					finally:
						if f != sys.stdout:
							f.close()
							
					badloops=0
				
				except KeyboardInterrupt:
					raise
				except:
					badloops+=1
					if badloops >= 5: raise #re-raise on 5th consecutive failure
	
		finally:
			lp.close()
		
	
	#monitor_humidity(filepath=None, sample_period=3.0)
	#monitor_humidity()
	
	if sys.platform=="darwin":
		if not "/sw/bin" in os.environ["PATH"]:
			os.environ["PATH"]=os.environ["PATH"]+":/sw/bin/"
		if not "DISPLAY" in os.environ:
			os.environ["DISPLAY"]=":0.0" #in case we want to use X and haven't set it up.
		
	from pyx import *
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
		
	
	def rt_test():
		badloops=0
		lp=LabPro_Mac_USB()
		
		try:	
			lp.wake()
			lp.reset()
			print lp.get_system_config() #just to vacuum the connection
			lp.setup_channel(chan=1, operation=1) #set up -10-10 volts readback
			lp.binary_mode(blocking_factor=4)
			 
			lp.setup_data_collection(samptime=0.0167, numpoints=-1, rectime=0)
			rtdata=[]
			starttime=time.time()
			for i in range(500):
				time.sleep(0.04)
				rtdata+=lp.get_data_binary_realtime()
			lp.stop()
			rtdata+=lp.get_data_binary_realtime()
			
			stoptime=time.time()
			
			rtdata=Numeric.array(rtdata,Numeric.Float)
			volts=lp.scale_binary_data(rtdata[:,0],(-10.,10.))
			timevals=Numeric.add.accumulate(rtdata[:,1]*0.0001)
			print stoptime-starttime, len(rtdata), Numeric.array_str(volts[:20],precision=4, suppress_small=1)
			print timevals[-1], Numeric.array_str(timevals[:20],precision=4, suppress_small=1)
	
			g =graphxy(width=20, y=graph.linaxis(title='volts (V) \\& current (mA)'), x=graph.linaxis(title='time (sec)', min=0) )
			g.texrunner.settex(lfs='foils17pt')
						
			g.plot(
				graph.data(data.data(zip(timevals,volts)), x=0, y=1),
				graph.line() )
				
			g.display('test')
	
	
		finally:
			lp.close()
			
	rt_test()
