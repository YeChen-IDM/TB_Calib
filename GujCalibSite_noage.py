import logging

from calibtool.study_sites.IncidenceCalibSite import IncidenceCalibSite
from TBTImeSeriesCalibSite import TBTimeSeriesCalibSite

logger = logging.getLogger(__name__)


class GujCalibSite( TBTimeSeriesCalibSite):

    reference_data = {

        # Taken from GTB 2015 estimates
        # produced by WHO last lookup September 28th 2017

        "Year Prevalence": [2011],
        "Active pulmonary prevalence": [ 0.00264 ],
        "Year Incidence": [2012, 2013, 2014, 2015, 2016],
        "Annual TB Incidence": [
            139, 148, 154, 164, 173],

    }

    def __init__(self):
        super(GujCalibSite, self).__init__('Gujarat')
