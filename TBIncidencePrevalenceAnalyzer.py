
# TODO: Generalize this style of CalibAnalyzer as much as possible
#       to minimize repeated code in e.g. PrevalenceByAgeAnalyzer
import logging
import os

import pandas as pd

from calibtool import LL_calculators
from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer
from tb.utils.TB_LL_calculators import LL_geometric_mean_normal as LL_geometric_mean_normal

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
        # incidence_variance_key = 'TB Incidence 95'
        #mortality_key = 'Annual TB Mortality'
        #mortality_variance_key = 'TB Mortality 95'
        if incidence_key in dictkeys:
            self.reference = pd.Series(site.reference_data.get(incidence_key) ).apply(lambda x: x/100000.0)
        # if incidence_variance_key in dictkeys:
        #     self.reference_variance = pd.Series(site.reference_data.get(incidence_variance_key)).apply(lambda x: (1/1.96 * (x/1e5))**2)
        # if mortality_key in dictkeys:
        #     self.reference_mortality = pd.Series(site.reference_data.get(mortality_key)).apply(lambda x: x / 100000.0)
        # if mortality_variance_key in dictkeys:
        #     self.reference_variance_mortality = pd.Series(site.reference_data.get(mortality_variance_key)).apply(
        #         lambda x: (1 / 1.96 * (x / 1e5)) ** 2)

        self.reference_years_prev = site.reference_data['Year Prevalence']

    def filter(self, sim_metadata):
        '''
        This analyzer only needs to analyze simulations for the site it is linked to.
        N.B. another instance of the same analyzer may exist with a different site
             and correspondingly different reference data.
        '''
        #return sim_metadata.get('__site__', False) == self.site.name
        return True

    def apply(self, parser):
        '''
        Extract data from output data
        '''
        data_yrs = self.site.reference_data['Year Prevalence']
        # data = pd.DataFrame.from_dict(parser.raw_data[self.filenames[0]])
        data = parser.raw_data[self.filenames[0]].copy()
        a = data['Year']
        data['Year'] = data['Year'].apply(lambda x: int(round(x)) + self.start_year)
        data = data[data['Year'].between(min(data_yrs),max(data_yrs))]

        data_reduced = data.groupby('Year').mean()
        data_reduced.reindex(['Year'])
        data_reduced = data_reduced.loc[data_yrs]

        channel_data = data_reduced
        index_number = parser.sim_data.get('__sample_index__')
        channel_data.sample = parser.sim_data.get('__sample_index__')                  #parser.sim_data.get('__sample_index__')
        channel_data.sim_id = parser.sim_id

        channel_data = pd.DataFrame({'sample': index_number,
                                    'sim_id': parser.sim_id,
                                    'TB_Incidence': data_reduced['Incidence'].tolist(),
                                    'TB_Prevalence': data_reduced['Active_Sx'].tolist(),
                                    'Year': data_reduced.index.values})
        channel_data.set_index(['sample','sim_id','Year'], drop= True,inplace= True)

        return channel_data

    def combine(self, parsers):
        '''
        Combine the simulation data into a single table for all analyzed simulations.
        '''
        selected = [p.selected_data[id(self)] for p in parsers.values() if id(self) in p.selected_data ]
        combined = pd.concat(selected, axis= 'index')
        #combined.groupby(['sample','sim_id'])['channels']['TB Incidence'].mean()
        #combined.reindex(['sample', 'sim_id'])
        print(combined)
        #self.data = combined.set_index(['sample','sim_id','Year'], drop= True,inplace= False)
        self.data = combined
        logger.debug(self.data)

    def compare(self, sample):
        '''
        Assess the result per sample, in this case the likelihood
        comparison between simulation and reference data.
        '''
        #print(sample)
        sample = sample.groupby(level = ['sample', 'Year']).mean()
        if len(sample.index.get_level_values('sample').unique()) != 1 :
            raise ValueError('There should only be one sample')
        shoot = sample.index.droplevel(level= 'sample').tolist()
        sample.reset_index(drop = True, inplace= True)

        LL_inc =  LL_geometric_mean_normal(self.reference,sample['TB_Incidence'],self.reference_variance)
        LL_mort = LL_geometric_mean_normal(self.reference_mortality,sample['TB_Mortality'],self.reference_variance_mortality)

        return (LL_inc + LL_mort)



    def finalize(self):
        '''
        Calculate the output result for each sample.
        '''
        print(self.data)
        self.result = self.data.groupby(level=['sample']).apply(self.compare)
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
