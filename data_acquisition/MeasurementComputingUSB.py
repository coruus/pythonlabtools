"MeasurementComputingUSB supports connections of  Measurement Computing, Inc.  USB devices"

_rcsid="$Id: MeasurementComputingUSB.py,v 1.2 2003-11-14 20:13:51 mendenhall Exp $"



class MeasurementComputingError(Exception):
	pass

import array
import Numeric

import time
import sys

import os
import threading
import traceback
import struct

from operator import isSequenceType

try:
	import fcntl #on platforms with fcntl, use it!
	
	class USB_libusb_mixin:
		"mixin class  to allow operation of MCC devices via USB port using libusb pipe server"
		
		server_executable_path=os.path.join(os.path.dirname(__file__),"MeasurementComputingServer")
		
		def __init__(self):
			self.usb_send, self.usb_recv, self.usb_err=os.popen3(self.server_executable_path+( " %d" % -self.device_index),'b',0)
			self.usb_timestamp=0
			self.__usb_read_leftovers=''
			try:
				fcntl.fcntl(self.usb_recv, fcntl.F_SETFL, os.O_NONBLOCK) #pipes must be nonblocking
				fcntl.fcntl(self.usb_err, fcntl.F_SETFL, os.O_NONBLOCK) #pipes must be nonblocking
				firstmsg=''
				badloops=0
				while badloops<5  and not firstmsg:
					time.sleep(0.5)
					try:
						firstmsg=self.usb_err.read()
					except:
						badloops+=1
					
				if badloops==5:
					raise MeasurementComputingError("cannot find or communicate with USB Server: "+str(self.server_executable_path))
				print firstmsg
				if firstmsg.find("Found device") < 0 :
					raise MeasurementComputingError("USB Server could not connect to a MCC device at index %d" % self.device_index)  
				
				self.idProduct=int(firstmsg[firstmsg.find("ID=")+3:],16) #which flavor of MCC device?
				self.__keep_running=1
				self.status_monitor=threading.Thread(target=self.read_status, name='USB MCC status monitor')
				self.status_monitor.start()
				self.__saved_realtime_fragment=''
			except:
				self.__keep_running=0
				self.close()
				raise
				
		def check_usb_status(self):
			if not self.__keep_running or not self.status_monitor.isAlive():
				raise MeasurementComputingError("MCC USB server has died...")
				
		def high_speed_serial(self):
			"not implemented on USB"
			return 
	
		def high_speed_setup(self):
			"not implemented onUSB"
			return 
	
		def set_port_params(self):
			"not implemented on USB"
			return 
	
		def read(self, maxlen=None):
			"read data from USB"
			self.check_usb_status()
			res=''
			
			db=self.__usb_read_leftovers
					
			while(not maxlen or (maxlen and len(res) < maxlen)):
				while len(db)<20: #8 bytes + 8 byte timestamp + 4 byte 0xffffffff flag		
					try:
						db+=self.usb_recv.read(20-len(db)) 
					except IOError:
						err=sys.exc_info()[1].args
						if err[0] in (11, 29, 35): #these errors are sometimes returned on a nonblocking empty read
							pass #just return empty data
						else:
							print  "USB server disconnected unexpectedly", err
							raise MeasurementComputingError("USB server disconnected unexpectedly", sys.exc_info()[1].args)
	
					if not maxlen and not db: break #no data at all, just fall out of this inner loop
	
				if not maxlen and not db: break #doing uncounted read, just take data until it quits coming
								
				flag, tv_sec, tv_usec=struct.unpack('LLL', db[:12])
				if flag != 0x00ffffff:
					raise MeasurementComputingError("Bad packet header from MCC device: " + ("%04x %08x %08x"%(flag, tv_sec, tv_usec)))
				self.data_timestamp=float(tv_sec)+float(tv_usec)*1e-6
				res+=db[12:20]
				db=db[20:]

			self.__usb_read_leftovers=db		
			
			return res
				
	
		def read_status(self):
			"monitor the pipe server's stderr in a thread"
			while(self.__keep_running):
				try:
					res=self.usb_err.read()
					if res: print >> sys.stderr, "MCC USB message: ", res
					if res.find("****EXITED****") >= 0: break #if this thread dies, LabPro has gone away somehow
				except IOError:
					if sys.exc_info()[1].args[0] in (11,29, 35): #this error is returned on a nonblocking empty read
						time.sleep(1) #so just wait and try again
					else:
						raise
			
		def write(self, data_array):
			self.check_usb_status()
			datastr="%d "*len(data_array) % tuple(data_array) + "\n"
			#print "writing...", data_array, datastr
			self.usb_send.write(datastr)
			#time.sleep(0.02)
						
		def close(self):
			self.__keep_running=0
			try:
				self.stop()
			except:
				pass
			try:
				self.usb_send.write("****QUIT****\n")
			except:
				pass
			time.sleep(2)
			self.usb_recv.close()
			self.usb_send.close()
			self.usb_err.close()	
		
	default_server_mixin=USB_libusb_mixin
	
except ImportError:
	
	class USB_Win32_libusb_mixin:
		"mixin class for RawLabPro to allow operation of LabPro via USB port on Win32 using pipe server"
		
		server_executable_path=os.path.join(os.path.dirname(__file__),"MeasurementComputingServer.exe")
		
		def __init__(self):
			self.usb_send, self.usb_recv, self.usb_err=os.popen3(self.server_executable_path+( " %d" % -self.device_index),'b',-1)
			self.usb_timestamp=0
			self.__usb_read_leftovers=''
			try:
				firstmsg=''
				badloops=0
				while badloops<5 and (not firstmsg or not firstmsg[-1] in '\r\n'):
					try:
						firstmsg+=self.usb_err.readline(1000)
					except:
						traceback.print_exc()
						badloops+=1
					
				if badloops==5:
					raise MeasurementComputingError("cannot find or communicate with USB Server: "+str(self.server_executable_path))
				if firstmsg.find("No LabPro Found") >=0:
					raise MeasurementComputingError("USB Server could not connect to an MCC device at index %d" % self.device_index)  
	
				self.__keep_running=1
				self.status_monitor=threading.Thread(target=self.read_status, name='USB MCC status monitor')
				self.status_monitor.start()
				self.__saved_realtime_fragment=''
			except:
				self.__keep_running=0
				self.close()
				raise
				
		def check_usb_status(self):
			if not self.__keep_running or not self.status_monitor.isAlive():
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
	
		def read(self, maxlen=None):
			"""read data from USB"""
			self.check_usb_status()
			res=''
			db=self.__usb_read_leftovers
					
			while( not maxlen or len(res) < maxlen):
				while len(db)<20: #8 bytes + 8 byte timestamp + 4 byte 0xffffffff flag
	
					#a pipe is a seek-and-tallable object, so we can see how much data is there this way			
					self.usb_recv.seek(0,2)
					count=self.usb_recv.tell()
					self.usb_recv.seek(0)
					
					if count: db+=self.usb_recv.read(count)
					if not db: break #no data at all, just fall out of this inner loop
	
				if not db: break #doing uncounted read, just take data until it quits coming
								
				flag, tv_sec, tv_usec=struct.unpack('LLL', db[:12])
				if flag != 0x00ffffff:
					raise LabProError("Bad packet header from MCC device: " + ("%04x %08x %08x"%(flag, tv_sec, tv_usec)))
				self.data_timestamp=float(tv_sec)+float(tv_usec)*1e-6
				res+=db[12:20]
				db=db[20:]
				
			self.__usb_read_leftovers=db		
			return res
				
	
		def read_status(self):
			"monitor the pipe server's stderr in a thread"
			while(self.__keep_running):
				res=self.usb_err.readline(1000)
				if res: print >> sys.stderr, "MCC USB message: ", res
				if res.find("****EXITED****") >= 0: break #if this thread dies, LabPro has gone away somehow
			
		def write(self, data_array):
			self.check_usb_status()
			datastr="%d "*len(data_array) % tuple(data_array) + "\n"
			#print "writing...", data_array, datastr
			self.usb_send.write(datastr)
			time.sleep(0.02) #give a little extra time for message passing
			
		def close(self):
			time.sleep(2)
			self.__keep_running=0
			try:
				self.stop()
			except:
				pass
			try:
				self.usb_send.write("****QUIT****\n")
			except:
				pass
			time.sleep(2)
			self.usb_recv.close()
			self.usb_send.close()
			self.usb_err.close()	
	
	default_server_mixin=USB_Win32_libusb_mixin

class MCC_Device(default_server_mixin):
	
	GAIN1_SE=0x08
	GAIN1_DIFF=0x00
	GAIN2_DIFF=0x10
	GAIN4_DIFF=0x20
	GAIN5_DIFF=0x30
	GAIN8_DIFF=0x40
	GAIN10_DIFF=0x50
	GAIN16_DIFF=0x60
	GAIN20_DIFF=0x70
	
	BASE_CLOCK_RATE=6000000 #correct for MiniLab1008
	TIMER_STEPS=( #prescaler and setup times for rates between given entry and next entry
			(100, 7, 0),
			(200, 6, 0),
			(400, 5, 0),
			(800, 4, 1),
			(1500, 3, 3),
			(3000, 2, 6),
			(6000, 1, 10),
			(8000, 0, 0) #no operation beyond this speed
	)
	
	DATA_SAMPLE_SIZE = 64
	DATA_PACKET_SIZE=96
	MAX_STORED_SAMPLES=4096
	
	AD_BURST_MODE=1
	AD_CONT_MODE=2
	AD_SINGLE_MODE=4
	AD_NO_CAL=8
	AD_EXTTRIGGER=16
		
	ad_gain_scale_dict = {
			GAIN1_SE :  (20.0/4096.0, 10.0, 1.0),
			GAIN1_DIFF: (40.0/4096.0, 20.0, 1.0),
			GAIN2_DIFF: (40.0/4096.0, 20.0, 2.0),
			GAIN4_DIFF: (40.0/4096.0, 20.0, 4.0),
			GAIN5_DIFF: (40.0/4096.0, 20.0, 5.0),
			GAIN8_DIFF: (40.0/4096.0, 20.0, 8.0),
			GAIN10_DIFF: (40.0/4096.0, 20.0, 10.0),
			GAIN16_DIFF: (40.0/4096.0, 20.0, 16.0),
			GAIN20_DIFF: (40.0/4096.0, 20.0, 20.0)
	}
		
	def counts_to_volts(self, gain, counts):
		if gain==self.GAIN1_SE:
			counts=counts*2
		else:
			counts = counts ^ 0x800 #sign bit is inverted, apparently, per MCC DLL
		scale1, offset, scale2 = self.ad_gain_scale_dict[gain]
		return ((counts*scale1)-offset)/scale2		
			
	def __init__(self, device_index=1):
		self.device_index=device_index
		default_server_mixin.__init__(self)
		self.scanning=0
		
	def blink_led(self):
		self.write((11,))
		time.sleep(2)
		
	def analog_output(self, channel, volts):
		counts=max(min(int(volts*4095.0/5.0 + 0.5), 4095),0)
		lowcounts=counts & 255
		highcounts = counts // 256
		self.write((8, channel, lowcounts, highcounts))
	
	def analog_input(self, channel, gain=GAIN1_DIFF):
		self.write((6,channel, gain))
		retdata=self.read(8)
		pair=ord(retdata[0]), ord(retdata[1])
		if retdata:
			retval=ord(retdata[0]) + ( ord(retdata[1]) << 4 )
			return self.counts_to_volts(gain, retval)
		else:
			return None

	def setup_gain_list(self, channels, gains):
		changain=[c & 0x07 | g for c,g in zip(channels, gains)]
		while len(channels) >=6:
			self.write([7,6]+changain[:6])
			changain=changain[6:]
		if changain:
			self.write([7,len(changain)]+changain)

	def compute_timer_vals(self, Rate):
		ts=self.TIMER_STEPS
		if Rate < ts[0][0] or Rate > ts[-1][0]:
			raise MeasurementComputingError, "sample rate not in range %d-%d" % (ts[0][0], ts[-1][0])
		
		#find appropriate range in list
		lower=ts[0]
		for entry in ts[1:]:
			if entry[0] > Rate: break
			lower=entry
		
		baserate, timerPre, setupTime=lower
		timerMult= 1 << (timerPre+1)
		timerVal= int((256 - (self.BASE_CLOCK_RATE/(Rate * timerMult))) + 0.5)
	
		actRate = (float(self.BASE_CLOCK_RATE)/((256-timerVal) * timerMult))
	
		return timerPre, timerVal, setupTime, actRate
		
	def setup_analog_scan(self, channels=(0,2,4,6), sweeps=-1, rate=100, gains=GAIN2_DIFF, exttrig=0):
	
		if self.scanning:
			raise MeasurementComputingError, "cannot start new scan without stopping old one"
			
		if not isSequenceType(channels):
			channels=(channels,)
		if not isSequenceType(gains):
			gains=(gains,)*len(channels)
			
		if len(gains) != len(channels):
			raise MeasurementComputingError, \
				"gain list not compatible with channel list" +str(channels)+":"+str(gains)
		
		self.setup_gain_list(channels, gains)

		timerPre, timerVal, setupTime, actRate=self.compute_timer_vals(rate*len(channels))
		
		if sweeps<=0: #continuous scan !
			blocking = self.DATA_SAMPLE_SIZE//len(channels)
			if len(channels)*blocking != self.DATA_SAMPLE_SIZE:
				raise MeasurementComputingError, \
					"continuous scan channel count must be submultiple of %d" % self.DATA_SAMPLE_SIZE
			tCount=self.DATA_SAMPLE_SIZE #each block is filled in continuous mode
			scanmode=self.AD_CONTINUOUS
			
		elif sweeps==1: #use AD_SINGLEEXEC mode for single sweep
			tCount=len(channels)
			scanmode=self.AD_SINGLE_MODE
			
		else:
			#round block count up to next block above requested samples
			blocks=(len(channels)*sweeps+self.DATA_SAMPLE_SIZE-1)//self.DATA_SAMPLE_SIZE
			tCount=blocks*self.DATA_SAMPLE_SIZE
			if tCount>self.MAX_STORED_SAMPLES:
				raise MeasurementComputingError, \
					"burst scan sample count must be <  %d" % self.MAX_STORED_SAMPLES
			scanmode=self.AD_BURST_MODE
				
		print locals()
		
		tChigh=tCount>>8
		tClow=tCount&255
		self.write((14, tClow, tChigh, timerVal+setupTime, timerPre, scanmode | exttrig))
	
	
if __name__=='__main__':
	
	mcc=MCC_Device()
	time.sleep(0.5)
	try:
		mcc.blink_led()
		
		if 0:
			for i in range(2):
				mcc.analog_output(0, 4.5*(i & 1))
				time.sleep(0.5)
			
			for i in range(200):
				print mcc.analog_input(0, gain=mcc.GAIN5_DIFF)
		elif 1:
			mcc.setup_analog_scan()
					
	finally:
		mcc.close()
		

	pass