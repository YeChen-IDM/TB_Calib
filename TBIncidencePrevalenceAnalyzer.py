
# TODO: Generalize this style of CalibAnalyzer as much as possible
#       to minimize repeated code in e.g. PrevalenceByAgeAnalyzer
import logging
import os

import pandas as pd

from calibtool import LL_calculators
from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer
from tb.utils.TB_LL_calculators import LL_geometric_mean_normal as LL_geometric_mean_normal
import math

logger = logging.getLogger(__name__)


class TBIncidencePrevalenceAnalyzer(BaseCalibrationAnalyzer):


    data_group_names = ['sample', 'sim_id', 'channels']

    def __init__(self, site, weight=1, name= 'TimeSeries' ):
        super(TBIncidencePrevalenceAnalyzer, self).__init__(site)
        self.name = name
        self.weight = weight
        self.filenames = [os.path.join('output','Report_TBHIV_ByAge.csv')]
        self.setup = {}
        self.start_year = 1700 # 2000
        self.set_site(site)

    def set_site(self, site):
        '''
        Get the reference data that this analyzer needs from the specified site.

        Get survey collection dates and subregions, if present, from the specified site.
        '''
        self.site = site
        dictkeys = site.reference_data.keys()
        incidence_key = 'Annual TB Incidence'
        incidence_variance_key = 'Annual TB Incidence_95'
        prevalence_key = "Active pulmonary prevalence"
        prevalence_variance_key = 'Active pulmonary prevalence_95'
        if incidence_key in dictkeys:
            self.reference_inci = pd.Series(site.reference_data.get(incidence_key) ).apply(lambda x: x/100000.0)
        if prevalence_key in dictkeys:
            self.reference_prev = pd.Series(site.reference_data.get(prevalence_key)).apply(lambda x: x / 100000.0)
        self.reference_variance_prev = self.reference_prev
        if incidence_variance_key in dictkeys:
            self.reference_variance_inci = pd.Series(site.reference_data.get(incidence_variance_key)).apply(lambda x: (1/1.96 * (x/1e5))**2)
        if prevalence_variance_key in dictkeys:
            self.reference_variance_prev = pd.Series(site.reference_data.get(prevalence_variance_key)).apply(
                lambda x: (1 / 1.96 * (x / 1e5)) ** 2)

        self.reference_years_prev = site.reference_data['Year Prev']
        self.reference_years_inci = site.reference_data['Year Incidence']


    def filter(self, sim_metadata):
        '''
        This analyzer only needs to analyze simulations for the site it is linked to.
        N.B. another instance of the same analyzer may exist with a different site
             and correspondingly different reference data.
        '''
        #return sim_metadata.get('__site__', False) == self.site.name
        return True

    def last(self, df):
        return df.ix[[-1]]

    def apply(self, parser):
        '''
        Extract data from output data
        '''
        data_yrs_prev = self.site.reference_data['Year Prev']
        data_yrs_inci = self.site.reference_data['Year Incidence']
        # data = pd.DataFrame.from_dict(parser.raw_data[self.filenames[0]])
        data = parser.raw_data[self.filenames[0]].copy()
        # a = data['Year']
        data['Year'] = data['Year'].apply(lambda x: math.floor(x) + self.start_year)
        data_prev = data[data['Year'].between(min(data_yrs_prev),max(data_yrs_prev))]
        data_inci = data[data['Year'].between(min(data_yrs_inci),max(data_yrs_inci))]

        agebin_map = {'LESS_1': '0-14', 'LESS_5': '0-14', 'LESS_10': '0-14','LESS_15': '0-14',
                                     'LESS_20': '15-19', 'LESS_25': '20-24', 'LESS_30': '25-29', 'LESS_35': '30-34',
                                     'LESS_40': '35-39', 'LESS_45': '40-44', 'LESS_50': '45-49', 'LESS_55': '50-54',
                                     'LESS_60': '55-59', 'LESS_65': '60-64', 'LESS_70': '65-69', 'LESS_75': '70-74', 'LESS_80': '75-79',
                                     'LESS_85': '80+', 'LESS_90': '80+', 'LESS_95': '80+', 'GREAT_95': '80+'}
        count = 0
        for x, _ in data_prev.iterrows():
            # drop the first half year's prevalence data
            if count < len(agebin_map):
                data_prev.drop(x, inplace=True)
            count += 1
            if count == len(agebin_map) * 2:
                count = 0
        for df in [data_prev, data_inci]:
            df.replace({'AgeBin': agebin_map}, inplace=True)
            # df.sort_values(['AgeBin','Year'], inplace=True)

        data_prev_reduced = data_prev.groupby(['AgeBin', 'Year'], as_index=True).sum()
        data_inci_reduced = data_inci.groupby(['AgeBin', 'Year'], as_index=True).sum()


        # data_reduced = data.groupby('Year').mean()
        # data_reduced.reindex(['Year'])
        # data_reduced = data_reduced.loc[data_yrs]
        #
        # channel_data = data_reduced
        # index_number = parser.sim_data.get('__sample_index__')
        # channel_data.sample = parser.sim_data.get('__sample_index__')                  #parser.sim_data.get('__sample_index__')
        # channel_data.sim_id = parser.sim_id
        index_number = parser.sim_data.get('__sample_index__')
        channel_data_prev = pd.DataFrame({'sample': index_number,
                                     'sim_id': parser.sim_id,
                                     'TB_Prevalence': ((data_prev_reduced['Active'] -
                                                        data_prev_reduced['PrevalentExtraPulmonary'])
                                                       /data_prev_reduced['Population']).tolist(),
                                     'agebin_year': data_prev_reduced.index.values})
        channel_data_inci = pd.DataFrame({'sample': index_number,
                                     'sim_id': parser.sim_id,
                                     'TB_Incidence': data_inci_reduced['Incidence'].tolist(),
                                     'agebin_year': data_inci_reduced.index.values})
        channel_data_prev.set_index(['sample','sim_id','agebin_year'], drop= True,inplace= True)
        channel_data_inci.set_index(['sample','sim_id','agebin_year'], drop= True,inplace= True)

        return [channel_data_prev, channel_data_inci]

    def combine(self, parsers):
        '''
        Combine the simulation data into a single table for all analyzed simulations.
        '''
        selected = [p.selected_data[id(self)] for p in parsers.values() if id(self) in p.selected_data ]
        combined_prev = pd.concat([selected_list[0] for selected_list in selected], axis= 'index')
        combined_inci = pd.concat([selected_list[1] for selected_list in selected], axis= 'index')
        #combined.groupby(['sample','sim_id'])['channels']['TB Incidence'].mean()
        #combined.reindex(['sample', 'sim_id'])
        #print(combined_prev, combined_inci)
        #self.data = combined.set_index(['sample','sim_id','Year'], drop= True,inplace= False)
        self.data = [combined_prev, combined_inci]
        logger.debug(self.data)

    def compare(self, sample, reference, reference_variance, column_name):
        '''
        Assess the result per sample, in this case the likelihood
        comparison between simulation and reference data.
        '''
        #print(sample)
        sample = sample.groupby(level=['sample', 'agebin_year']).mean()
        if len(sample.index.get_level_values('sample').unique()) != 1:
            raise ValueError('There should only be one sample')
        #shoot = sample.index.droplevel(level='sample').tolist()
        sample.reset_index(drop=True, inplace=True)

        LL = LL_geometric_mean_normal(reference, sample[column_name], reference_variance)

        return LL

    def finalize(self):
        '''
        Calculate the output result for each sample.
        '''
        print(self.data)
        self.result = self.data[0].groupby(level=['sample']).apply(
            self.compare, reference=self.reference_prev, reference_variance=self.reference_variance_prev,
            column_name='TB_Prevalence') \
                      + self.data[1].groupby(level=['sample']).apply(
            self.compare,  reference=self.reference_inci, reference_variance=self.reference_variance_inci,
            column_name='TB_Incidence')
        logger.debug(self.result)

    def cache(self):
        '''
        Return a cache of the minimal data required for plotting sample comparisons
        to reference comparisons.
        '''
        #cache = self.data.copy()

        #sample_dicts = []
        #for idx, df in cache.groupby(level='sample', sort=True) :
        #    d = { 'region' : self.regions,
        #           self.y : [sdf[self.y].values.tolist() for jdx, sdf in df.groupby(level='region') ] }
        #    sample_dicts.append(d)

        #logger.debug(sample_dicts)

        #return {'sims': sample_dicts, 'reference': self.reference, 'axis_names': ['region', self.y]}
        return []

    def uid(self):
        ''' A unique identifier of site-name and analyzer-name. '''
        return '_'.join([self.site.name, self.name])

    def plot_comparison(cls, fig, data, **kwargs):
        """
        Plot data onto figure according to logic in derived classes
        """
        pass
