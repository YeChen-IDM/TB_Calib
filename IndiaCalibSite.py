import logging

from calibtool.study_sites.IncidenceCalibSite import IncidenceCalibSite
from TBTImeSeriesCalibSite import TBTimeSeriesCalibSite

logger = logging.getLogger(__name__)


class IndiaCalibSite(TBTimeSeriesCalibSite):

    reference_data = {
        "Years": [2000, 2001, 2002, 2003, 2004, 2005, 2006,
                  2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015],
        "Annual TB Incidence": [
            585, 666, 746, 820, 883,
            932, 963, 977, 977, 967,
            948, 922, 892, 860, 834,
            834],
        "Annual TB Mortality": [
            158, 196, 179, 202, 200,
            263, 251, 249, 226, 210,
            202, 196, 206, 191, 178,
            179],
        "TB Incidence 95": [
            229, 260, 291, 320, 345,
            364, 376, 381, 279, 255,
            254, 223, 273, 268, 260,
            327],
        "TB Mortality 95": [
          90, 95, 123, 131, 149,
          39, 39, 41, 37, 46,
          57, 73, 79, 93, 109,
          104]
    }

    def __init__(self):
        super(IndiaCalibSite, self).__init__('Gujarat')
