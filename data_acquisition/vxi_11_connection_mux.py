"""connection_multiplexer allows dynamic opening and closing of vxi_11 connections to (probably low-priority) devices, 
so a huge number of static connections isn't needed."""
_rcsid = "$Id$"

import vxi_11
import time
import threading


class connection_mux:

  "multiplex vxi_11 connections of all types through a single handler.  <name> must be globally unique."

  # keys used in dictionaries
  # a dictionary key for the max number of active connections on a given host
  max_connections = "_max"
  # a dictionary key for the current number of active connections on a given
  # host
  current_connections = "_current"
  # a dictionary key for the details about a specific vxi_11 host
  host_info = "_hostinfo"
  # a dictionary key for the list of active supplies on a given host
  active_connections = "_active"
  host_lock = "_hostlock"

  default_max_active_connections = 3

  def __init__(self):
    self.devices = {}
    self.hosts = {}

  def add_host(self, hostname, hostinfo={}, max_connections=None):
    "add a host"
    if max_connections is None:
      max_connections = self.default_max_active_connections
    self.hosts[hostname] = {
        self.host_info: hostinfo, self.max_connections: max_connections, self.
        active_connections: {},
        self.host_lock: threading.Condition()}

  def add_device_and_init(
          self, name, hostname, device_class, address, channel_parameters):
    """add a device by class and initialize it. 
    channel_parameters contains whatever info is needed by the device to tell one channel from another
    """
    physical_key = name  # this should be unique for a given physical device on a given host
    if not physical_key in self.hosts[hostname]:
      device = device_class(
          name, self.hosts[hostname][self.host_info],
          address)  # never seen this one before, initialize it
    self.add_device(name, device, channel_parameters, physical_key)

  def add_device(self, name, device, channel_parameters=None,
                 physical_key=None):
    "add an already initialized device, and insert appropriate replacements for lock_connection and unlock_connection"
    if not physical_key:
      # this should be unique for a given physical device on a given host
      physical_key = name
    hostname = device.host  # device should have this field available
    if not physical_key in self.hosts[hostname]:
      self.hosts[hostname][physical_key] = device
      # we may not be connected, because of previous mux activity
      if device.connected:
        # don't hang onto the connection now, we will reconnect when we need it
        device.disconnect()
      device.setup_mux(self, name)

    self.devices[name] = (hostname, physical_key, channel_parameters)

  def __getitem__(self, name):
    hostname, physical, params = self.devices[name]
    return self.hosts[hostname][physical]

  def __setitem__(self, name, device):
    self.add_device(name, device)

  def lock_connection(self, name):
    "dynamically allocate and lock  a connection to the named device"
    hostname, physical, params = self.devices[name]
    host = self.hosts[hostname]
    lock = host[self.host_lock]
    lock.acquire()
    try:
      supply = host[physical]
      conns = host[self.active_connections]
      # connection is not already established, must negotiate
      if not physical in conns:
        # we need to allocate a connection
        while len(conns) >= host[self.max_connections]:
          # and we need to throw one out first, so find oldest unused key
          # note that a key points to [lockflag, time], and lockflag=1 always
          # sorts higher than 0, so unlocked connections sort first
          k = conns.keys()
          kmin = k[0]
          tmin = conns[kmin]
          # print "****  "+str(conns)
          for i in k[1:]:
            # comparing lists, so first element (locked) has priority, then
            # last-use time
            if conns[i] < tmin:
              kmin = i
              tmin = conns[i]

          # best device is locked, so none are available to deallocate!
          if tmin[0]:
            lock.wait()  # get in line for a slot when unlock lets one go
          else:
            # mark as in-transit-closing as an error flag for later debugging
            tmin[0] = 2
            # terminate this connection (very unlikely to fail)
            host[kmin].disconnect()
            del conns[kmin]  # and delist it if it closed correctly

        supply.reconnect()  # try to reconnect
        conns[physical] = [1, time.time(), 0]  # take a new queue slot
      else:  # we are still in queue, so grab the lock
        conns[physical][0] = 1  # marked locked
        conns[physical][2] = 0  # reset to default priority
        # note: we haven't disturbed the activity timer here
        # this way, if someone set a high priority somewhere else, it won't be
        # forgotten
    finally:
      lock.release()

  def release_all(self):
    "really close all connections we are still holding.  Why would we want to do this?"
    for name in self.devices:
      hostname, physical, params = self.devices[name]
      host = self.hosts[hostname]
      device = host[physical]
      # OK to get this outside of lock, since it is mutable and  so is always
      # up to date
      conns = host[self.active_connections]
      lock = host[self.host_lock]
      lock.acquire()
      if physical in conns:  # it is really connected
        if conns[physical][0]:  # oops, device is still locked
          device.log_error("Attempt to release locked device from mux")
        else:
          try:
            device.disconnect()
          except:
            device.log_error("Cannot release device from mux!?")
        del conns[physical]

  def adjust_priority(self, name, priority=0):
    "set priority so next real unlock will see it.  Increases trump decreases"
    if not priority:
      return  # default priority doesn't change anything else
    hostname, physical, params = self.devices[name]
    host = self.hosts[hostname]
    lock = host[self.host_lock]
    lock.acquire()
    # paranoid... nothing in here can fail, but use safe programming style
    try:
      conn = host[self.active_connections][physical]
      conn2 = conn[2]
      conn[2] = max(priority, conn2 - 10000 * (conn2 == 0))
      # use max(priority, conn2) unless conn2==0 (defaulted), in which case
      # priority always wins
    finally:
      lock.release()

  def unlock_connection(self, name, priority=0):
    "dynamically mark unlockable  a connection to the named device, and set (biased) last-used time to time.time()+priority"
    hostname, physical, params = self.devices[name]
    host = self.hosts[hostname]
    lock = host[self.host_lock]
    lock.acquire()
    try:
      supply = host[physical]
      conn = host[self.active_connections][physical]
      conn[0] = 0  # reset the lock flag
      # adjust apparent last-used time, but only carefully,
      # so non-default settings, especially positive ones, are 'sticky'
      expire, conn2 = conn[1:3]
      now = time.time()
      delta = expire - now  # current remaining priority
      # the logic:  new clock = now + highest of non-defaulted (priorities and remaining time on timer)
      # where remaining time on timer defaults if timer has expired, and
      # priorities default if = 0
      actpri = max(
          priority - 10000 * (priority == 0),
          conn2 - 10000 * (conn2 == 0),
          delta * (delta > 0) - 10000 * (delta <= 0))
      # print "** %12s %4.0f %4.0f %5.1f %.1f" % (name, priority, conn2, delta,
      # actpri)
      # new expiration is biased only if actpri isn't defaulted to -10000
      newtime = now + actpri * (actpri != -10000)
      conn[1] = newtime
    finally:
      lock.notify()  # in case someone is waiting for a slot
      lock.release()

  def delete_physical_device(self, name):
    "delete all channels attached to the same device as <name> from the mux"
    hostname, physical, params = self.devices[name]
    host = self.hosts[hostname]
    lock = host[self.host_lock]
    lock.acquire()
    try:
      conns = host[self.active_connections]
      if physical in conns:
        # device better be unlocked!
        assert not conns[physical][
            0], "Attempt to delete device from mux while locked!"
        # mark as in-transit-closing as an error flag for later debugging
        conns[physical][0] = 2
        host[physical].close()
        del conns[physical]
      host[physical].setup_mux(None, None)  # release device from this mux
      del host[physical]
      # remove all names attached to this physical device on this host
      for dev in self.devices:
        if self.devices[dev][:2] == [hostname, physical]:
          del self.devices[dev]
    finally:
      lock.release()
