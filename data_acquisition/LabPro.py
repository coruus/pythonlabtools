"""LabPro supports communications with the Vernier Instruments (www.vernier.com) LabPro Module
over a serial line"""

rcsid="$Id: LabPro.py,v 1.1 2003-05-23 17:40:28 mendenhall Exp $"

import time
import Numeric
import os
import sys
import math
import termios
import array

import struct
_bigendian=(struct.unpack('H','\00\01')[0] == 1)

class LabProError(Exception):
	"general class for LabPro Errors"
	pass

class LabProTimeout(LabProError):
	"class for LabPro transient errors, which might be recoverable"
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
	_highspeedports.append(port_name)
		
class LabPro:
	commands={'reset':0, 'channel_setup':1, 'data_collection_setup':3, 
			'conversion_equation_setup':4, 'data_control':5, 'system_setup':6, 
			'request_system_status':7, 'request_channel_status':8, 'request_channel_data':9,
			'digital_data_capture':12, 'baudrate':105 , 'request_channel_setup':115,
			'analog_output':401, 'led':1998}

	def __init__(self, port_name, highspeed=0):
		
		init_port_memory()
		self.setup_serial(port_name)
		self.__highspeed=0
		self.__saved_realtime_frag=''
		if port_name in _highspeedports:
			self.high_speed_serial()
		elif highspeed:
			self.high_speed_setup()
			remember_highspeed_port(port_name)
		self.wake()
		self.wake()
		
	def set_port_params(self, baud=termios.B38400):
		port=self.serial_read_port
		attrs=termios.tcgetattr(port)
		attrs[4]=attrs[5]=baud #set 38.4kbaud
		attrs[2] = attrs[2] | termios.CLOCAL #ignore connection indicators
		cc=attrs[6]; cc[termios.VMIN]=0; cc[termios.VTIME]=0 #non-blocking reads
		termios.tcsetattr(port, termios.TCSADRAIN, attrs)
	
	def setup_serial(self, port_name):
		self.serial_write_port=port=open(port_name,"wb", 0)
		self.serial_read_port=open(port_name,"rb", 16384)
		self.set_port_params() #setup default port
	
	def high_speed_serial(self):
		"use this if you now the LabPro is already running at high speed"
		if self.__highspeed:	return  #already did this
		try:
			self.set_port_params(baud=termios.B115200) #see if this is allowed!
		except:
			return
		#if this isn't an error, try to do it for real!	
		self.set_port_params(baud=termios.B115200) #should work now, since it worked before
		time.sleep(0.1)
		self.__highspeed=1

	def high_speed_setup(self):
		"use this to command the LabPro up to high speed, then switch the serial up. Do not use if it has already been commanded!"
		self.wake()
		self.command('baudrate',115)
		time.sleep(0.5)
		self.serial_read_port.read(10)
		self.high_speed_serial()
		
	def send_string(self, s='s'):
		self.serial_write_port.write(s+'\r')
		termios.tcdrain(self.serial_write_port)
		time.sleep(0.05)
		
	def command(self, name, *parms):	
		s='s{%d' %  self.commands[name]
		for p in parms:
			s=s+','+str(p)
		s=s+'}'
		self.send_string(s)
	
	def read_ascii_response(self):
		str=''
		empties=0
		while(empties < 5 and str[-3:]!='}\r\n'):
			time.sleep(.1)
			newdata=self.serial_read_port.read(1000)
			str+=newdata
			if newdata:
				empties=0
			else:
				empties+=1	
		if empties: #last result must have gotten data, so empties should be zero
			raise LabProTimeout('timeout getting ascii data, current result: '+repr(str))		
		goodstart=str.find('{')
		assert goodstart>=0, 'bad ascii data: '+repr(str)
		return map(float,str[goodstart+1:-3].split(','))


	def reset(self):
		self.send_string('s')
		self.command('reset')
		time.sleep(0.25) #extra-long sleep


	def reset_highspeed_serial(self):
		self.__highspeed=0

	def close(self):
		self.serial_read_port.close()
		self.serial_write_port.close()
		
	def wake(self):
		self.send_string('s')

	def get_system_config(self):
		self.command('request_system_status')
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
		self.__save_realtime_frag='' #just for safety
		self.command('data_collection_setup', 
				samptime, numpoints, trigtype, trigchan, trigthresh, prestore, 0, 
				rectime, filter, fastmode)
		
	def setup_channel(self, chan=0, operation=1, postproc=0, equation=0):
		self.command('channel_setup', chan, operation, postproc, equation)

	def wait_for_data_done(self, flasher=0):
		while(1):
			if flasher:
				self.flash_led('yellow', 0.25)
			state=self.get_system_config()
			if state['error']:
				raise LabProError(state['error'])
			syst=state['system state'] & 7
			if syst==1:
				raise LabProError('Waiting for data when none being collected')
			if syst==4:
				break #all done
		return state['sample count']

	def labpro_get_data_ascii(self, chan=0):
		labpro_command('data_control',chan,0,0,0,1)
		labpro_send_string('g')
		return self.read_ascii_response()

	def labpro_get_data_ascii_realtime():
		s=self.__saved_realtime_frag
		while(s.find('}\r\n') < 0):
			time.sleep(0.1)
			s=s+serial_read_port.read()
			if not s:
				return [] #no data if string is still blank
		sl=s.split('}\r\n')
		if s[-3:] != '}\r\n': #last number was incomplete
			self.__saved_realtime_frag=sl[-1]
		else:
			self.__saved_realtime_frag=''
		return [map(float, ptstr[1:].split(',')) for ptstr in sl[:-1]]

	def binary_mode(self):
		self.command('conversion_equation_setup', 0,-1)


	def parse_binary(self, data, channels, chunklen):
		#print repr(data)
		chk=reduce(lambda x,y: x^y, array.array('B',data))
		assert chk ==0 or chk==0xff, "Bad checksum on string: "+('%02x '%chk)+repr(data)
		format='>'+channels*'H'+'L' #data+time+check
		l=[struct.unpack(format,data[i:i+chunklen-1]) for i in range(0,len(data), chunklen)]
		return l
		
	def get_data_binary_realtime(self, channels=1):
		s=self.__saved_realtime_frag
		chunklen=2*channels+5 #data+timestamp(4 bytes)+chkbyte
		while(len(s)<chunklen):
			time.sleep(0.1)
			s=s+self.serial_read_port.read()
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
		scaled_min, scaled_max=scaled_range	
		return Numeric.array(idata,Numeric.Float)*((scaled_max-scaled_min)/65536.0)+scaled_min
	
	def get_data_binary(self, chan=0, points=None, scaled_range=None):
		if points is None:
			points=self.get_system_config()['dataend']
		self.command('data_control',chan,0,0,0,1)
		self.send_string('g')
		chunklen=2*points+1
		s=''
		empties=0
		while(empties<5 and len(s)<chunklen):
			time.sleep(0.1)
			newdata=self.serial_read_port.read(chunklen-len(s))
			s+=newdata
			if newdata:
				empties=0
			else:
				empties+=1
		
		if empties:
			raise LabProTimeout('timeout getting binary data, got %d bytes, expected %d, stuff looks like: %s' %(len(s), chunklen, repr(s[:100])))		
			
		chk=reduce(lambda x,y: x^y, array.array('B',s))
		if not (chk ==0 or chk==0xff):
			raise LabProError( "Bad checksum on string: "+('%02x '%chk)+repr(s[:100]))
		data=array.array('H')
		data.fromstring(s[:-1])
		if not _bigendian:
			data.byteswap()
		if scaled_range is None:
			return data #return raw integers	
		else:
			return self.scale_binary_data(data, scaled_range) 
	
	def stop(self):
		self.command('system_setup',0)

	def get_current_channel_data(self, channel=1):
		self.command('request_channel_data', channel, 0)
		self.send_string('g')
		returnself.read_ascii_response()[0]

	def get_channel_setup(self, channel=1):
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
		if waveform=='dc':
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
			offset=int(0.5+(amplitude*0.5-center)/0.0006)
	
		self.command('analog_output', mode, ampl, offset, int(1000.0*period_sec))

	def analog_output_off():
		self.command('analog_output',0,0,0,0)

	ledcolors={'red':1, 'yellow':2, 'green':3}
	
	def set_led(self, color='red', state=1):
		led=self.ledcolors[color.lower()]
		self.command('led', led, state)
	
	def flash_led(self, color='green', flashtime=0.5):
		self.set_led(color,1)
		time.sleep(flashtime)
		self.set_led(color,0)
		
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
