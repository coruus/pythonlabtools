/* serve up USB data from a Vernier LabPro device attached via USB using libusb on MacOSX, Linux or *BSD */

static char rcsid[]="RCSID $Id: LabProlibusbServer.c,v 1.10 2003-11-06 19:56:47 mendenhall Exp $";

/* 
requires libusb (from www.sourceforge.net) installed 
to compile on a Mac under OSX:
cc -o LabProUSBServer -framework IOKit -framework CoreFoundation -lusb LabProlibusbServer.c

or to compile it using libusb statically, (so it can be installed on machines without libusb)
cc -o LabProUSBServer -framework IOKit -framework CoreFoundation /usr/local/lib/libusb.a LabProlibusbServer.c

to compile under Linux or *BSD:
cc -o LabProUSBServer -lusb LabProlibusbServer.c

that should produce a working binary.
*/

#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <pthread.h>
#include <signal.h>

#include <usb.h>

int keep_running=1;
usb_dev_handle *global_intf; /* global reference to device, for cleanup */
int use_time_stamps=0;

void handle_signal(int what)
{
	int err;
	keep_running=0;
	if (global_intf) {
		usb_reset(global_intf); /* terminate eternal read operation */
		sleep(1);
		usb_reset(global_intf); /* terminate eternal read operation */
	}
	fprintf(stderr,"Got signal\n");
	
}

int pass_input(usb_dev_handle *udev)
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
			count = usb_bulk_write(udev, USB_ENDPOINT_OUT | 2, buf, currbufpos, 1000);
			if (count < 0 || count != currbufpos)
			{
				fprintf(stderr, "write error: count=%d,  %s\n", count, usb_strerror());
				break;
			}
			currbufpos=0;
		}
	}
	keep_running=0;
	return 0;
}


int pass_output(usb_dev_handle *udev)
{
	int err, count;
	const unsigned int retbufsize=64; /* LabPro always transfers 64 bytes blocks */
	struct { int blockflag; struct timeval tv; char inBuf[retbufsize];} datastruct;
	struct timezone tz;
	time_t start_time, stop_time;
	
	datastruct.blockflag=0x00ffffff; /* make it easy to find timestamps in data */
	
	while(keep_running) {
		count=0;
		start_time=time(NULL);
		count = usb_bulk_read(udev, USB_ENDPOINT_IN | 2 , datastruct.inBuf, retbufsize, 1000000);
		if (keep_running && count != retbufsize) {
			stop_time=time(NULL);
			if(stop_time-start_time < 1000) {
				/* timeouts are 1000 seconds (1000000 milliseconds), so if we fail after this long, it's
					probably a timeout */
				fprintf(stderr, "read error: %s\n", usb_strerror());
				break;
			} else continue;
		} 
		if(keep_running) {
			if(use_time_stamps) {
				gettimeofday(&datastruct.tv, &tz);
				err=write(fileno(stdout), (void *)&datastruct, sizeof(datastruct));
			} else {
				err=write(fileno(stdout), &datastruct.inBuf, retbufsize);
			}
			fflush(stdout);
			if (err<0) break;
		}
	}
	keep_running=0;
	return 0;	
}

void dealWithDevice(struct usb_device *dev, usb_dev_handle *udev)
{
	int err=1,i;
	pthread_t input_thread, output_thread;
	void *thread_retval;
	
	/* sometime other processes may be probing the LabPro just when we try to claim it, so try a few times */
	for(i=0; i<3 || err; i++) {	
		err=usb_claim_interface(udev, 0);
		if(err) sleep(1);
	}
	if (err) {
		fprintf(stderr, "error: %s\n", usb_strerror());
		return;
	}
	
	if (dev->config) {
		/* this is what should be done, but the LabPro has no descriptors, so we will set the value to 1 if dev->config is NULL */
		usb_set_configuration(udev, dev->config[0].bConfigurationValue); /* configure interface */
	} else {
		usb_set_configuration(udev, 1); /* configure interface */
	} 
	
	err=pthread_create(&input_thread, 0, (void *)pass_input, udev);
	if(!err) err=pthread_create(&output_thread, 0, (void *)pass_output, udev);
	
	if(!err) {
		err=pthread_join(input_thread, &thread_retval);
		usb_reset(udev); /* terminate eternal read operation */
		sleep(1);
		usb_reset(udev); /* terminate eternal read operation */
		err=pthread_join(output_thread, &thread_retval);
	}
	usb_reset(udev);
	usb_release_interface(udev,0);
	
}


int main (int argc, const char * argv[])
{
    int			idVendor = 0x8f7;
    int			idProduct = 1;
    int USBIndex, matchcount;
	int i;
	usb_dev_handle *udev=0;
	struct usb_bus *bus;
	struct usb_device *dev, *matchdev;
	
	/* if one argument is provided, it should be an index as to _which_ USB LabPro is to be opened 
		providing a negative index enables time stamping as is used when this is run as a robot */
	if (argc==2) {
		USBIndex=atoi(argv[1]);
		if (abs(USBIndex) < 1 || abs(USBIndex) > 255) {
			fprintf(stderr,"Bad USB index argument provided... should be 1<=index<=255 or negative to enable binary time stamps, got: %s\n", argv[1]);
			fprintf(stderr,"****EXITED****\n");
			return 1;
		}
		if (USBIndex < 0) {
			USBIndex=-USBIndex;
			use_time_stamps=1;
		}
		USBIndex -=1;
	} else USBIndex=0;
	
	setbuf(stdout, 0);
	setbuf(stderr,0);

	usb_init();
	usb_find_busses();
	usb_find_devices();
    
	signal(SIGHUP, handle_signal);
	signal(SIGINT, handle_signal);
	signal(SIGQUIT, handle_signal);
	signal(SIGTERM, handle_signal);
	signal(SIGPIPE, handle_signal);
	
	matchcount=-1;
	for (bus = usb_busses; bus && matchcount<USBIndex; bus = bus->next) {
		for (dev = bus->devices; dev && matchcount<USBIndex; dev = dev->next) {			
			if(dev->descriptor.idVendor==idVendor && dev->descriptor.idProduct==idProduct) {
				matchcount++; matchdev=dev; }
		}
	}
	

	if(matchcount==USBIndex) {
		udev = usb_open(matchdev);
		if(udev) {
			fprintf(stderr, "Found device %p\n", (void*)udev);
			fflush(0);
			global_intf=udev;
			usb_set_debug(0);
			dealWithDevice(matchdev, udev);
			global_intf=0; /* don't need resets any more */
			usb_reset(udev);
			usb_close(udev);
		} else {
			fprintf(stderr, "Found but couldn't open device %d... probably already open\n", USBIndex+1);
		}
    } else fprintf(stderr,"No LabPro Found at index %d\n", USBIndex+1);
	    
	fprintf(stderr,"****EXITED****\n");
    return 0;
}
