#!/usr/bin/env python

#dddbtag = "dddb-20130929"
#dddbtag = "dddb-20130929-1"
#dddbtag = "dddb-20150928"

dddbtag = "dddb-20130929"

from Configurables import LHCbApp
LHCbApp().DDDBtag = dddbtag


# =============================================================================
import PartProp.PartPropAlg
import PartProp.Service
from   GaudiPython.Bindings import AppMgr
# =============================================================================

gaudi = AppMgr()

gaudi.initialize()

pps   = gaudi.ppSvc()

bsprop = pps.find('B_s0')
bslprop = pps.find('B_s0L')
bshprop = pps.find('B_s0H')

print '\nBs lifetime        = %.3f ps\n' % (bsprop.lifeTime() * 1e3)
print '\nBsL lifetime        = %.3f ps\n' % (bslprop.lifeTime() * 1e3)
print '\nBsH lifetime        = %.3f ps\n' % (bshprop.lifeTime() * 1e3)
