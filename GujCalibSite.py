import logging

from calibtool.study_sites.IncidenceCalibSite import IncidenceCalibSite
from TBTImeSeriesCalibSite import TBTimeSeriesCalibSite

logger = logging.getLogger(__name__)


class GujCalibSite( TBTimeSeriesCalibSite):

    reference_data = {

        # Age ranges 0-14, 15-19, 20-24, 25-29, 30-34, 35-39, 40-44, ...., 75-79, 80+
        # Incidence matrix: each row is yearly data for that age range
        "Year Prev": [2011],
        "Active pulmonary prevalence": [0.00053, 0.00284, 0.00292, 0.00227, 0.00301, 0.00359, 0.00272, 0.00610, 0.00469, 0.00462, 0.00503, 0.00663, 0.00530, 0.00796],
        "Year Incidence": [2012, 2013, 2014, 2015, 2016],
        "Annual TB Incidence": [
            20, 20, 21, 22, 21,
            99, 113, 116, 124, 134,
            186, 202, 207, 223, 225,
            189, 199, 207, 224, 230,
            194, 209, 216, 221, 226,
            213, 226, 228, 236, 251,
            222, 233, 237, 252, 257,
            244, 253, 269, 289, 312,
            281, 303, 320, 336, 367,
            220, 236, 245, 273, 298,
            322, 335, 362, 407, 461,
            225, 242, 262, 286, 326,
            198, 217, 238, 238, 291,
            122, 128, 148, 162, 209,
            127, 152, 144, 186, 187,
        ],
    }

    def __init__(self):
        super(GujCalibSite, self).__init__('Gujarat')
