/* serve up USB data from a Vernier LabPro device attached to a macintosh via USB */

static char _rcsid[]="$Id: LabProUSBMacServer.c,v 1.2 2003-06-04 17:37:29 mendenhall Exp $";

/* to compile on a Mac under OSX:
cc -o LabProUSBMacServer -framework IOKit -framework CoreFoundation LabProUSBMacServer.c

that should produce a working binary.
*/

#include <mach/mach.h>
#include <CoreFoundation/CFNumber.h>
#include <IOKit/IOKitLib.h>
#include <IOKit/IOCFPlugIn.h>
#include <IOKit/usb/IOUSBLib.h>

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>

mach_port_t 	masterPort = 0;				// requires <mach/mach.h>

int keep_running=1;

int pass_input(IOUSBInterfaceInterface182 **intf)
{
	fd_set inpipeinfo;
	struct timeval timeout;
	int hasdata, count, err;
	char buf[1024];
	static int currbufpos=0;

	while(keep_running) {
		timeout.tv_sec=0;
		timeout.tv_usec=10;
	
		FD_ZERO(&inpipeinfo);
		FD_SET(fileno(stdin), &inpipeinfo);
		count=select(1, &inpipeinfo, 0, 0, &timeout);
		hasdata=FD_ISSET(fileno(stdin), &inpipeinfo);
		if(!(count && hasdata)) continue; /* select says no data on stdin */
		count=read(fileno(stdin), &buf[currbufpos], sizeof(buf)-10-currbufpos);	
		if (count <=0) continue; /* somehow, even though select saw data, this didn't ?! */
		currbufpos+=count;
		buf[currbufpos]=0;
		if(buf[currbufpos-1]=='\r' || buf[currbufpos-1]=='\n') { /* detect line ending */
			if (strncmp(buf,"****QUIT****",12)==0) break;
			err = (*intf)->WritePipeTO(intf, 2, buf, currbufpos, 1000, 2000);
			currbufpos=0;
			if (err)
			{
				fprintf(stderr, "transferData: WritePipeAsyncFailed, err = %08x\n", err);
				break;
			}
		}
	}
	keep_running=0;
	return 0;
}

int pass_output(IOUSBInterfaceInterface182 **intf)
{
	int err;
	char inBuf[128];
	UInt32 retbufsize=64; /* LabPro always transfers 64 bytes blocks */

	while(keep_running) {
		err = (*intf)->ReadPipeTO(intf, 1, inBuf, &retbufsize, 1000, 2000);
		
		if (err==kIOReturnTimeout || err==(int)0xe0004051) continue; /* no data, just go on */
		else if (err) {
			fprintf(stderr, "transferData: ReadPipeASync, err = %08x\n", err);
			break;
		}
		err=write(fileno(stdout), inBuf, 64);
		fflush(stdout);
		if (err<0) break;
	}
	keep_running=0;
	return 0;	
}

void transferData(IOUSBInterfaceInterface182 **intf)
{
	int err=0;
	pthread_t input_thread, output_thread;
	void *thread_retval;
	
	/*
	char startup[]="s\ns{0}\ns{7}\n";

	err = (*intf)->WritePipeTO(intf, outPipeRef, startup, strlen(startup), 1000, 2000);
	if (err)
	{
		fprintf(stderr, "transferData: WritePipeAsyncFailed, err = %08x\n", err);
		return;
	}
	
	*/
	
	err=pthread_create(&input_thread, 0, (void *)pass_input, intf);
	if(!err) err=pthread_create(&output_thread, 0, (void *)pass_output, intf);
	
	if(!err) {
		err=pthread_join(input_thread, &thread_retval);
		err=pthread_join(output_thread, &thread_retval);
	}
}

void dealWithInterface(io_service_t usbInterfaceRef)
{
    IOReturn				err;
    IOCFPlugInInterface 		**iodev;		// requires <IOKit/IOCFPlugIn.h>
    IOUSBInterfaceInterface182 		**intf;
    SInt32 				score;
    UInt8				numPipes;


    err = IOCreatePlugInInterfaceForService(usbInterfaceRef, kIOUSBInterfaceUserClientTypeID, kIOCFPlugInInterfaceID, &iodev, &score);
    if (err || !iodev)
    {
		fprintf(stderr, "dealWithInterface: unable to create plugin. ret = %08x, iodev = %p\n", err, iodev);
	return;
    }
    err = (*iodev)->QueryInterface(iodev, CFUUIDGetUUIDBytes(kIOUSBInterfaceInterfaceID182), (LPVOID)&intf);
    (*iodev)->Release(iodev);				// done with this
    if (err || !intf)
    {
		fprintf(stderr, "dealWithInterface: unable to create a device interface. ret = %08x, intf = %p\n", err, intf);
	return;
    }
    err = (*intf)->USBInterfaceOpen(intf);
    if (err)
    {
		fprintf(stderr, "dealWithInterface: unable to open interface. ret = %08x\n", err);
	return;
    }
	transferData(intf);
		
    err = (*intf)->USBInterfaceClose(intf);
    if (err)
    {
		fprintf(stderr, "dealWithInterface: unable to close interface. ret = %08x\n", err);
	return;
    }
    err = (*intf)->Release(intf);
    if (err)
    {
		fprintf(stderr, "dealWithInterface: unable to release interface. ret = %08x\n", err);
	return;
    }
}

void dealWithDevice(io_service_t usbDeviceRef)
{
    IOReturn				err;
    IOCFPlugInInterface 		**iodev;		// requires <IOKit/IOCFPlugIn.h>
    IOUSBDeviceInterface 		**dev;
    SInt32 				score;
    UInt8				numConf;
    IOUSBConfigurationDescriptorPtr	confDesc;
    IOUSBFindInterfaceRequest		interfaceRequest;
    io_iterator_t			iterator;
    io_service_t			usbInterfaceRef;
    
    err = IOCreatePlugInInterfaceForService(usbDeviceRef, kIOUSBDeviceUserClientTypeID, kIOCFPlugInInterfaceID, &iodev, &score);
    err = (*iodev)->QueryInterface(iodev, CFUUIDGetUUIDBytes(kIOUSBDeviceInterfaceID), (LPVOID)&dev);
    (*iodev)->Release(iodev);				// done with this
    err = (*dev)->USBDeviceOpen(dev);
    if (err)
    {
		fprintf(stderr, "dealWithDevice: unable to open device. ret = %08x\n", err);
	return;
    }
    err = (*dev)->GetConfigurationDescriptorPtr(dev, 0, &confDesc);			// get the first config desc (index 0)
    if (err)
    {
		fprintf(stderr, "dealWithDevice:unable to get config descriptor for index 0\n");
        (*dev)->USBDeviceClose(dev);
        (*dev)->Release(dev);
	return;
    }
    err = (*dev)->SetConfiguration(dev, confDesc->bConfigurationValue);
    if (err)
    {
		fprintf(stderr, "dealWithDevice: unable to set the configuration\n");
        (*dev)->USBDeviceClose(dev);
        (*dev)->Release(dev);
	return;
    }
    
    interfaceRequest.bInterfaceClass = kIOUSBFindInterfaceDontCare;		// requested class
    interfaceRequest.bInterfaceSubClass = kIOUSBFindInterfaceDontCare;		// requested subclass
    interfaceRequest.bInterfaceProtocol = kIOUSBFindInterfaceDontCare;		// requested protocol
    interfaceRequest.bAlternateSetting = kIOUSBFindInterfaceDontCare;		// requested alt setting
    
    err = (*dev)->CreateInterfaceIterator(dev, &interfaceRequest, &iterator);    
    usbInterfaceRef = IOIteratorNext(iterator); /* get first and only interface */
	dealWithInterface(usbInterfaceRef);
	IOObjectRelease(usbInterfaceRef);				// no longer need this reference    
    IOObjectRelease(iterator);
    iterator = 0;

    err = (*dev)->USBDeviceClose(dev);
    err = (*dev)->Release(dev);
}

int main (int argc, const char * argv[])
{
    kern_return_t		err;
    CFMutableDictionaryRef 	matchingDictionary = 0;		// requires <IOKit/IOKitLib.h>
    SInt32			idVendor = 0x8f7;
    SInt32			idProduct = 1;
    CFNumberRef			numberRef;
    io_iterator_t 		iterator = 0;
    io_service_t		usbDeviceRef;
    
	setbuf(stdout, 0);
	setbuf(stderr,0);
	
    IOMasterPort(MACH_PORT_NULL, &masterPort);				
    matchingDictionary = IOServiceMatching(kIOUSBDeviceClassName);	// requires <IOKit/usb/IOUSBLib.h>
    numberRef = CFNumberCreate(kCFAllocatorDefault, kCFNumberSInt32Type, &idVendor);
    CFDictionaryAddValue(matchingDictionary, CFSTR(kUSBVendorID), numberRef);
    CFRelease(numberRef);
    numberRef = 0;
    numberRef = CFNumberCreate(kCFAllocatorDefault, kCFNumberSInt32Type, &idProduct);
    CFDictionaryAddValue(matchingDictionary, CFSTR(kUSBProductID), numberRef);
    CFRelease(numberRef);
    numberRef = 0;
    
    err = IOServiceGetMatchingServices(masterPort, matchingDictionary, &iterator);
    matchingDictionary = 0;			// this was consumed by the above call
    
	usbDeviceRef = IOIteratorNext(iterator); /* get first instance */
	if(usbDeviceRef) {
		fprintf(stderr, "Found device %p\n", (void*)usbDeviceRef);
		fflush(0);
		dealWithDevice(usbDeviceRef);
		IOObjectRelease(usbDeviceRef);			// no longer need this reference
    } else fprintf(stderr,"No LabPro Found\n");
	
    IOObjectRelease(iterator);
    iterator = 0;
    
    mach_port_deallocate(mach_task_self(), masterPort);
	fprintf(stderr,"****EXITED****\n");
    return 0;
}
