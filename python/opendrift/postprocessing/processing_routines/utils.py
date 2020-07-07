# $Id: utils.py 118 2009-01-20 17:54:48Z iant $
#
# Shared utilities.

import math

def formatDegrees(x):
    degrees = int(x+0.00000000001)
    minutes = ((x+0.00000000001) - degrees)*60.0
    ret = u"%d\u00b0" % degrees
    ret += "%.0f'" % minutes
    return ret

def formatDegreesE(x, pos):
    return formatDegrees(x) + 'E'

def formatDegreesNone(x, pos):
    return formatDegrees(x)

def formatDegreesW(x, pos):
    return formatDegrees(-1.0*x) + 'W'

def formatDegreesN(x, pos):
    return formatDegrees(x) + 'N'

def formatDegreesS(x, pos):
    return formatDegrees(-1.0*x) + 'S'

def formatDegreesNull(x, pos):
    return ''

def formatLog(x, pos):
    # Dodgy botch.
    if x == 1.0:
        return '1'
    else:
        dp = int(-math.log10(x)*1.00001)
        return '0.' + '0'*(dp-1) + '1'

def formatNull(x, pos):
    return ''
