import logging
from collections import OrderedDict
from abc import ABCMeta
import os

from calibtool.CalibSite import CalibSite
from calibtool.study_sites.site_setup_functions import \
    config_setup_fn, summary_report_fn, add_treatment_fn, site_input_eir_fn
from calibtool.analyzers.ChannelBySeasonAgeDensityCohortAnalyzer import ChannelBySeasonAgeDensityCohortAnalyzer
from calibtool.analyzers.Helpers import season_channel_age_density_json_to_pandas,\
    season_channel_age_density_csv_to_pandas
# from TBPrev_IncidenceMortalityAnalyzer import TBPrev_IncidenceMortalityAnalyzer
# from TBIncidenceMortalityAnalyzer import TBIncidenceMortalityAnalyzer
from TBIncidencePrevalenceAnalyzer import TBIncidencePrevalenceAnalyzer

#from tb_csv_timeseries import TBTimeseriesAnalyzer

logger = logging.getLogger(__name__)


class TBTimeSeriesCalibSite(CalibSite):
    """
    An abstract class that implements the simulation setup for TB time series
    - South Africa Country model
    - Nigeria Country Model
    """

    __metaclass__ = ABCMeta

    metadata = {
        'Data_Years': []
    }

    def get_setup_functions(self):
        return []

    def get_reference_data(self, reference_type):
        site_ref_type = 'tb_time_series'

        if reference_type is not site_ref_type:
            raise Exception("%s does not support %s reference_type, only %s.",
                            self.__class__.__name__, reference_type, site_ref_type)

    def get_analyzers(self):
        #return [TBPrev_IncidenceMortalityAnalyzer(self,weight=1, name="Gujarat")]
        # return [TBIncidenceMortalityAnalyzer(self, weight=1, name="Gujarat")]
        return [TBIncidencePrevalenceAnalyzer(self, weight=1, name="Gujarat")]