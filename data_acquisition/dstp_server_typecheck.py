from nati_dstp_basics import *
import tagged_data

data=tagged_data.tagged_data_system()
s=NATI_DSTP_DataServer()

data.define_field("test1", s, "funny_field", [1,2,3,'hello'], writable=1)

print data["test1"]

data.test1=[4,5,6,"goodbye"]

print data.test1

data.test1='foobar'


	 
