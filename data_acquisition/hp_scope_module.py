"module importing a scope for global access"
#$Id: hp_scope_module.py,v 1.2 2003-04-11 17:51:02 mendenhall Exp $

from vxi_11_scopes import hp54542

class my_hp54542(hp54542):
	pass
	
		
if 0:
	scope=my_hp54542(host="***REMOVED***", portmap_proxy_host="127.0.0.1", portmap_proxy_port=1111,
			 device="gpib0,7",  timeout=4000, device_name="hp54542",
			raise_on_err=1)
else:
	scope=my_hp54542(host="***REMOVED***", 
			 device="gpib0,7",  timeout=5000, device_name="hp54542",
			raise_on_err=1)

