/* serve up USB data from a Vernier LabPro device attached via USB using libusb on MacOSX, Linux or *BSD */

static char rcsid[]="RCSID $Id: LabProlibusbServer.c,v 1.1 2003-07-11 19:22:26 mendenhall Exp $";

/* 
requires libusb (from www.sourceforge.net) installed 
to compile on a Mac under OSX:
cc -o LabProUSBServer -framework IOKit -framework CoreFoundation -lusb LabProlibusbServer.c

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

void handle_signal(int what)
{
	int err;
	keep_running=0;
	if (global_intf) {
		err = usb_resetep(global_intf, USB_ENDPOINT_IN | 2); /* terminate read operation */
		err = usb_clear_halt(global_intf, USB_ENDPOINT_IN | 2); /* terminate read operation */
		err = usb_resetep(global_intf, USB_ENDPOINT_OUT | 2); /* terminate write operation */
		err = usb_clear_halt(global_intf, USB_ENDPOINT_OUT | 2); /* terminate write operation */
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
	char inBuf[128];
	unsigned int retbufsize=64; /* LabPro always transfers 64 bytes blocks */

	while(keep_running) {
		count=0;
		count = usb_bulk_read(udev, USB_ENDPOINT_IN | 2 , inBuf, retbufsize, 0);
		
		if (keep_running && count != retbufsize) {
			fprintf(stderr, "read error: %s\n", usb_strerror());
			break;
		} 
		if(keep_running) {
			err=write(fileno(stdout), inBuf, 64);
			fflush(stdout);
			if (err<0) break;
		}
	}
	keep_running=0;
	return 0;	
}

void dealWithDevice(usb_dev_handle *udev)
{
	int err;
	pthread_t input_thread, output_thread;
	void *thread_retval;
	
	err=usb_claim_interface(udev, 0);
    if (err)
    {
		fprintf(stderr, "error: %s\n", usb_strerror());
	return;
    }
	
	err=pthread_create(&input_thread, 0, (void *)pass_input, udev);
	if(!err) err=pthread_create(&output_thread, 0, (void *)pass_output, udev);
	
	if(!err) {
		err=pthread_join(input_thread, &thread_retval);
		err = usb_resetep(udev,1); /* terminate eternal read operation */
		err=pthread_join(output_thread, &thread_retval);
	}
	usb_reset(udev);
	usb_release_interface(udev,0);
	
}


int main (int argc, const char * argv[])
{
    kern_return_t		err;
    int			idVendor = 0x8f7;
    int			idProduct = 1;
    int USBIndex, matchcount;
	int i;
	usb_dev_handle *udev=0;
	struct usb_bus *bus;
	struct usb_device *dev, *matchdev;
	
	/* if one argument is provided, it should be an index as to _which_ USB LabPro is to be opened */
	if (argc==2) {
		USBIndex=atoi(argv[1]);
		if (USBIndex < 1 || USBIndex > 255) {
			fprintf(stderr,"Bad USB index argument provided... should be 1<=index<=255, got: %s\n", argv[1]);
			fprintf(stderr,"****EXITED****\n");
			return 1;
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
		fprintf(stderr, "Found device %p\n", (void*)udev);
		fflush(0);
		global_intf=udev;
		usb_set_debug(0);
		dealWithDevice(udev);
		usb_reset(udev);
		usb_close(udev);
    } else fprintf(stderr,"No LabPro Found at index %d\n", USBIndex+1);
	    
	fprintf(stderr,"****EXITED****\n");
    return 0;
}
