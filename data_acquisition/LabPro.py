"""LabPro supports communications with the Vernier Instruments (www.vernier.com) LabPro Module
over a serial line"""

_rcsid="$Id: LabPro.py,v 1.14 2003-06-05 19:03:04 mendenhall Exp $"

import time
import Numeric
import os
import sys
import math
import array
import types
import operator

import struct
_bigendian=(struct.unpack('H','\00\01')[0] == 1)

class LabProError(Exception):
	"general class for LabPro Errors"
	pass

class LabProTransientError(LabProError):
	"class for LabPro transient errors, which might be recoverable"
	pass
		
class LabProTimeout(LabProTransientError):
	"class for timeouts"
	pass

class LabProDataError(LabProTransientError):
	"class for badly formatted data, may be recoverable on a retry"
	pass
	
def init_port_memory():
	global _highspeedports
	#_highspeedports is a persistent memory of which serial ports have had their LabPros set to high speed
	try:
		a=_highspeedports
	except:
		_highspeedports=[]

def unset_highspeed_port(port_name):
	"remove a LabPro from the list of devices known to be operating at high speed"
	try:
		_highspeedports.remove(port_name)
	except:
		pass

def remember_highspeed_port(port_name):
	"mark a  LabPro as already being set for high-speed operation"
	init_port_memory()
	if port_name not in _highspeedports:
		_highspeedports.append(port_name)
		
class RawLabPro:
	"support for Vernier Software LabPro data acquisition device.  RawLabPro class has no I/O and needs a mixin for the serial port support"
	
	commands={'reset':0, 'channel_setup':1, 'data_collection_setup':3, 
			'conversion_equation_setup':4, 'data_control':5, 'system_setup':6, 
			'request_system_status':7, 'request_channel_status':8, 'request_channel_data':9,
			'digital_data_capture':12, 'baudrate':105 , 'request_channel_setup':115,
			'analog_output':401, 'led':1998, 'sound': 1999, 'dig_out':2001,
			'power_control':102}

	def __init__(self, port_name, highspeed=0):
		init_port_memory()
		self.setup_serial(port_name)
		self.__saved_realtime_frag=''
		if port_name in _highspeedports:
			self.high_speed_serial()
		elif highspeed:
			self.high_speed_setup()
			remember_highspeed_port(port_name)
		self.wake()
		self.wake()
		for i in range(3): #try to wake up comunications by getting one good interchange
			try:
				config=self.get_system_config()
			except:
				if i==2: raise #last try failed, allow exception out
			else:
				break #transfer worked, all done
						
	def high_speed_setup(self):
		"use this to command the LabPro up to high speed, then switch the serial up. Do not use if it has already been commanded!"
		self.write('s\rs\rs\rs\r')
		self.command('baudrate',115)
		self.command('baudrate',115)
		self.command('baudrate',115)
		time.sleep(0.5)
		self.read()
		self.high_speed_serial()
		self.write('s\rs\rs\rs\r')
		self.command('baudrate',115)
		self.command('baudrate',115)
		self.command('baudrate',115)
		self.read()
		
	def send_string(self, s='s', delay=0.05):
		"send a string with appropriate end-of-line, and wait the specified time"
		self.write(s+'\r')
		time.sleep(delay)
		
	def command(self, name, *parms, **kwargs):
		"send a command based on either a command number of a symbolic name, with specified parameters and optional delay=xxx"
		if type(name) is not types.IntType: #if command is not an int, look it up
			name=self.commands[name]
			
		s='s{%d' %  name
		for p in parms:
			s=s+','+str(p)
		s=s+'}'
		try:
			delay=kwargs['delay']
		except:
			delay=0.05
		self.send_string(s, delay)
	
	def read_ascii_response(self):
		"read an ascii numeric list reponse that looks like {1.23e4, 5.67e8,  ...}"
		str=''
		empties=0
		while(empties < 5 and str[-3:]!='}\r\n'):
			time.sleep(.1)
			newdata=self.read()
			str+=newdata
			if newdata:
				empties=0
			else:
				empties+=1	
		if empties: #last result must have gotten data, so empties should be zero
			raise LabProTimeout('timeout getting ascii data, current result: '+repr(str))		
		goodstart=str.find('{')
		if goodstart<0:
			raise LabProDataError('bad ascii data: '+repr(str))
		return map(float,str[goodstart+1:-3].split(','))


	def reset(self):
		"reset LabPro to powerupd state, except baud rate" 
		self.send_string('s')
		self.command('reset', delay=0.25)

	def wake(self):
		"wake up LabPro"
		self.send_string('s', delay=0.01)

	def get_system_config(self):
		"return LabPro system status as a nice dictionary"
		self.command('request_system_status', delay=0.1)
		l=self.read_ascii_response()
		if  l[3]!=8888:
			raise LabProError("Bad system status from LabPro: "+str(l))
		dict={}
		
		dict['version']="%.5f"%l[0] #make string version
		for i, name in ((1,'error'), (2,'battery'),(5,'trigger condition'), (6,'channel function'),
				(7,'channel post'), (8,'channel filter'),(9,'sample count'),
				(10,'record time'),(12,'piezo flag'),(13,'system state'),
				(14,'datastart'), (15,'dataend'), (16,'system id')):
			dict[name]=int(l[i])
	 	
		for i, name in ((4,'sample time'),(11,'temperature')):
			dict[name]=float(l[i])
		
		return dict 

	def setup_data_collection(self, samptime=0.1, numpoints=100, trigtype=0, trigchan=0,
			trigthresh=0.0, prestore=0, rectime=0, filter=0, fastmode=0):
		"set up and start LabPro data acquisition"
		self.__saved_realtime_frag='' #just for safety
		self.command('data_collection_setup', 
				samptime, numpoints, trigtype, trigchan, trigthresh, prestore, 0, 
				rectime, filter, fastmode, delay=0.1)
		
	def setup_channel(self, chan=0, operation=1, postproc=0, equation=0):
		"configure a LabPro channel"
		self.command('channel_setup', chan, operation, postproc, equation)

	def wait_for_data_done(self, flasher=0):
		"wait until LabPro indicates data acquisition is complete"
		misses=0
		while(misses<5):
			if flasher:
				self.flash_led('yellow', 0.25)
			try:
				state=self.get_system_config()
			except:
				misses+=1
				continue
			
			misses=0
			if state['error']:
				raise LabProError(state['error'])
			syst=state['system state'] & 7
			if syst==1:
				raise LabProError('Waiting for data when none being collected')
			if syst==4:
				break #all done
	
		if misses:
			raise LabProError("5 consecutives retry failures in wait_for_data_done()... bailing out")
		
		return state['sample count']

	def get_data_ascii(self, chan=0):
		"tell LabPro to return the data from one channel"
		self.command('data_control',chan,0,0,0,1)
		self.send_string('g')
		return self.read_ascii_response()

	def get_data_ascii_realtime(self):
		"read the most recently collected realtime data from the LabPro and return as a list"
		s=self.__saved_realtime_frag
		while(s.find('}\r\n') < 0):
			time.sleep(0.1)
			s=s+self.read()
			if not s:
				return [] #no data if string is still blank
		sl=s.split('}\r\n')
		if s[-3:] != '}\r\n': #last number was incomplete
			self.__saved_realtime_frag=sl[-1]
		else:
			self.__saved_realtime_frag=''
		return [map(float, ptstr[1:].split(',')) for ptstr in sl[:-1]]

	def binary_mode(self):
		"tell LabPro that next data transfers will be in binary mode"
		self.command('conversion_equation_setup', 0,-1)

	def parse_binary(self, data, channels, chunklen):
		"convert LabPro binary string to integer list"
		#print repr(data)
		chk=reduce(operator.xor, array.array('B',data))
		if (chk!=0 and chk!=0xff):
			raise LabProDataError("Bad checksum on string: "+('%02x '%chk)+repr(data))
		format='>'+channels*'H'+'L' #data+time+check
		l=[struct.unpack(format,data[i:i+chunklen-1]) for i in range(0,len(data), chunklen)]
		return l
		
	def get_data_binary_realtime(self, channels=1):
		"read the most recently collected realtime data from the LabPro and return as a list, using binary transfers and no scaling"
		s=self.__saved_realtime_frag
		chunklen=2*channels+5 #data+timestamp(4 bytes)+chkbyte
		while(len(s)<chunklen):
			time.sleep(0.2)
			s=s+self.read()
			if not s:
				return [] #no data if string is still blank
		leftovers=len(s) % chunklen
		if leftovers:
			self.__saved_realtime_frag=s[-leftovers:]
			s=s[:-leftovers]
		else:
			self.__saved_realtime_frag=''
	
		return self.parse_binary(s, channels, chunklen)

	def scale_binary_data(self, idata, scaled_range):
		"scale data from LabPro binary transfer to specified full-scale min and max and return in Numeric (NumPy) array"
		scaled_min, scaled_max=scaled_range	
		return Numeric.array(idata,Numeric.Float)*((scaled_max-scaled_min)/65536.0)+scaled_min
	
	def get_data_binary(self, chan=0, points=None, scaled_range=None):
		"get specified channel data, either as raw integers if scaled_range=None, or as floats scaled to specified range in Numeric array"
		if points is None:
			points=self.get_system_config()['dataend']
		self.command('data_control',chan,0,0,0,1)
		self.send_string('g')
		chunklen=2*points+1
		s=''
		empties=0
		while(empties<5 and len(s)<chunklen):
			time.sleep(0.1)
			newdata=self.read(chunklen-len(s))
			s+=newdata
			if newdata:
				empties=0
			else:
				empties+=1
		
		if empties:
			raise LabProTimeout('timeout getting binary data, got %d bytes, expected %d, stuff looks like: %s' %(len(s), chunklen, repr(s[:100])))		
			
		chk=reduce(operator.xor, array.array('B',s))
		if not (chk ==0 or chk==0xff):
			raise LabProDataError( "Bad checksum on string: "+('%02x '%chk)+repr(s[:100]))
		data=array.array('H')
		data.fromstring(s[:-1])
		if not _bigendian:
			data.byteswap()
		if scaled_range is None:
			return data #return raw integers	
		else:
			return self.scale_binary_data(data, scaled_range) 
	
	def stop(self):
		"abort LabPro data acquisition"
		self.command('system_setup',0, delay=0.25)
		
	def get_current_channel_data(self, channel=1):
		"read current value of a channel"
		self.command('request_channel_data', channel, 0)
		self.send_string('g')
		return self.read_ascii_response()[0]

	def get_channel_setup(self, channel=1):
		"get configuration of a channel in a convenient dictionary"
		self.command('request_channel_setup', channel)
		l= self.read_ascii_response()
		dict={}
		for i, name in ((0,'cbl_2_sigfigs'), (1,'sigfigs'),(6,'typical_nsamples'), (7,'typical_command'),
				(8,'equation_index'), (13,'cal_pages'),(14,'current_cal_page')):
			dict[name]=int(l[i])
	 	
		for i, name in ((2,'ymin'),(3,'ymax'),(4,'yscale'), (5,'sample_rate'), 
				(9,'warmup_time'), (10,'K0'), (11,'K1'), (12,'K2')):
			dict[name]=float(l[i])
		
		return dict 

	def analog_output(self, waveform='dc', center=0.0, amplitude=1.0, period_sec=1.0):
		"output the specified waveform, amplitude and center value on channel 4.  Beware many funny rules from tech manual about limitations"
		if waveform=='off':
			self.command('analog_output',0,0,0,0)
			return			
		elif waveform=='dc':
			mode=1
			if center>=0.0:
				offset=0
				ampl=int(0.5+center/0.0024)
			else:
				offset=int(0.5-center/0.0012)
				ampl=0
		elif waveform[:3]=='sin':
			mode=6
			ampl=(4,3,2,2,1,1,0,0,0)[int(2.0*amplitude+0.01)]
			offset=int(0.5+((4,2,1,0.5,0.25)[ampl]*0.75-center)*1024)	
		else:
			mode={'rampup':2, 'rampdown':3, 'triangle':4, 'square':5}[waveform] 
			ampl=int(0.5+amplitude/0.0012)
			offset=int(0.5+(amplitude-center)/0.0012)
	
		self.command('analog_output', mode, ampl, offset, int(1000.0*period_sec), delay=0.1)

	def analog_output_off(self):
		"turn off analog output"
		self.command('analog_output',0,0,0,0)

	#note an LED has a number of different names
	ledcolors={'red':1, 'RED':1, 'r':1, 'R':1, 1:1,
			'yellow':2, 'YELLOW':2,'y':2,'Y':2,2:2,
			'green':3, 'GREEN':3,'g':3,'G':3,3:3 }
	
	def set_led(self, color='red', state=1):
		"set specified LED on (state=1) or off (state=0)"
		led=self.ledcolors[color]
		self.command('led', led, state, delay=0.01)
	
	def flash_led(self, color='green', flashtime=0.5):
		"flash specified LED on and then off"
		self.set_led(color,1)
		time.sleep(flashtime)
		self.set_led(color,0)
	
	def pattern_led(self,pattern, delay=0.0):
		"output a list of LED states"
		for ops in pattern:
			self.set_led(ops[0],ops[1])
			time.sleep(delay)
			
	red_green_ripple=(
		('red',1),('yellow',1),('red',0),('green',1),('yellow',0),('green',0))
	
	def sound(self, frequency=500.0, duration=0.1):
		"make a beep"
		self.command('sound', int(0.5+duration*1e4), int(0.5+4e4/frequency))

	def dig_out(self, *points):
		"output listed values to digital out port"
		if type(points) is type.IntType: #convert int to tuple
			points=(points, )
		self.command('dig_out',*points)

try:
	import termios
	class termios_mixin:
		"mixin class for RawLabPro to support Unix, MacOSX and Linux termios serial control"
					
		def close(self):
			self.serial_read_port.close()
			if self.serial_write_port != self.serial_read_port:
				self.serial_write_port.close()
			
		def set_port_params(self, baud=termios.B38400):
			port=self.serial_read_port
			attrs=termios.tcgetattr(port)
			attrs[4]=attrs[5]=baud #set 38.4kbaud
			attrs[2] = attrs[2] | termios.CLOCAL #ignore connection indicators
			cc=attrs[6]; cc[termios.VMIN]=0; cc[termios.VTIME]=0 #non-blocking reads
			termios.tcsetattr(port, termios.TCSADRAIN, attrs)
		
		def setup_serial(self, port_name):
			self.serial_write_port=port=open(port_name,"r+" , 0)
			self.serial_read_port=self.serial_write_port	
			self.set_port_params() #setup default port
		
		def high_speed_serial(self):
			"use this if you know the LabPro is already running at high speed"
			self.set_port_params(baud=termios.B115200) #should work now, since it worked before
			time.sleep(0.1)
	
		def write(self, data):
			"override this if communication is not over normal serial"
			self.serial_write_port.write(data)
		
		def read(self, maxlen=None):
			"override this as for write()"
			if maxlen is None:
				return self.serial_read_port.read()
			else:
				return self.serial_read_port.read(maxlen)

	class LabPro(termios_mixin, RawLabPro):
		"default LabPro uses system native serial, on MacOSX & Linux, this is termios"
		pass
	
except: #no termios serial support
	#insert check for other system architectures here, such as Windows or MacOS classic
	pass

try:
	import vxi_11
	class e5810_serial_mixin:
		"mixin class for RawLabPro to allow operation of LabPro via remote network connection over Agilent E5810 interface"
		def setup_serial(self,port_name=None):
			self.vxi11_host=vxi_11.vxi_11_connection(host=port_name, device="COM1", raise_on_err=1, timeout=5000, 
					device_name="Serial port on E5810")
			self.vxi11_host.lock() #we are exclusive owners by default!
		

		def high_speed_serial(self):
			"not implemented on vxi-11"
			return 
	
		def high_speed_setup(self):
			"not implemented on vxi-11"
			return 
	
		def set_port_params(self):
			"not implemented on vxi-11"
			return 
	
		def read(self, maxlen=None):
			return self.vxi11_host.read(count=maxlen)[-1] #last element of read is actual data
	
		def write(self, data):
			self.vxi11_host.write(data)
	
		def close(self):
			self.vxi11_host.disconnect()

	class LabPro_e5810(e5810_serial_mixin, RawLabPro):
		pass
		
except: #apparently vxi-11 e5810 serial support can't be found, just ignore
	pass

if __name__=='__main__':
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
	
	def try_labpro(binary=1):
		
		lp=LabPro(highspeed=1, port_name='/dev/cu.USA28X213P1.1')
		
		try:
			print lp.get_system_config()
			lp.reset()
			lp.setup_channel(chan=1)
			lp.setup_channel(chan=2, operation=2)
			if 0:
				lp.setup_data_collection(samptime=0.001, numpoints=200)
				lp.wait_for_data_done()
				print map(lambda x: "%.4f" % x, lp.get_data(chan=1)[:100] )
				print map(lambda x: "%.4f" % x, lp.get_data(chan=2)[:100] )
			lp.setup_channel(chan=2, operation=0)
			if 0:		
				lp.setup_data_collection(samptime=0.02, numpoints=-1, rectime=0)
				data=[]
				for i in range(10):
					time.sleep(1)
					data+=lp.get_data_ascii_realtime()
				lp.stop()
				data+=lp.get_data_ascii_realtime()
				
				data=Numeric.array(data,Numeric.Float32)
				print Numeric.array_str(data[:,0],precision=4, suppress_small=1)
			if 0:		
				lp.binary_mode()
				lp.setup_data_collection(samptime=0.002, numpoints=-1, rectime=2)
				data=[]
				for i in range(10):
					time.sleep(0.1)
					data+=lp.get_data_binary_realtime(channels=1)
				labpro_stop()
				data+=lp.get_data_binary_realtime(channels=1)
				data=Numeric.array(data,Numeric.Int32)
				print data[:,0]
			if 1:
				lp.reset()
				lp.analog_output(waveform='sine', center=-0.25, amplitude=4.0, period_sec=0.05)
				lp.setup_channel(chan=4) #let channel read back its own data
				lp.setup_channel(chan=2, operation=22) #current monitor (hidden command)
				
				if binary:
					lp.binary_mode()
				
				
				deltat=0.0002; npoints=1000
				lp.setup_data_collection(samptime=deltat, numpoints=npoints)
				lp.wait_for_data_done()
				if  binary:
					current=lp.get_data_binary(chan=2, points=npoints, scaled_range=(-0.1, 0.1) )
					voltage=lp.get_data_binary(chan=4, points=npoints, scaled_range=(-10.0, 10.0) )
				else:
					current=Numeric.array(lp.get_data_ascii(chan=2))
					voltage=Numeric.array(lp.get_data_ascii(chan=4))
				
				print Numeric.array_str(current, precision=4, suppress_small=1)
				print Numeric.array_str(voltage, precision=4, suppress_small=1)
				
					
		finally:
			lp.stop()
			try:
				print lp.get_system_config()
			except:
				print 'cannot get final status'
			lp.reset()
			lp.close()
				
	try_labpro()
