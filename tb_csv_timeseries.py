import logging
import os

import pandas as pd

from dtk.utils.analyzers import default_select_fn, default_group_fn, default_filter_fn
from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from dtk.utils.analyzers.plot import plot_by_channel
from calibtool.analyzers.BaseCalibrationAnalyzer import BaseCalibrationAnalyzer

logger = logging.getLogger(__name__)


def default_plot_fn(df, ax, df2):
    grouped = df.groupby(level=['group'], axis=1)
    grouped2 = df.groupby(level=['group'], axis=1)
    m2 = grouped2.mean()
    a = m2.columns
    b= df2.columns
    m = grouped.mean()
    m.plot(ax=ax, legend=True)

class TBTimeseriesAnalyzer(BaseCalibrationAnalyzer):

    plot_name = 'ChannelPlots'
    data_group_names = ['group', 'sim_id', 'channel']
    ordered_levels = ['channel', 'group', 'sim_id']
    output_file = 'timeseries.csv'

    def __init__(self,site=[],name=[],weight=1, filename=os.path.join('output', 'TBOutput.csv'), filter_function=default_filter_fn,
                 select_function=default_select_fn, group_function=default_group_fn, plot_function=default_plot_fn,
                 channels=['Statistical Population' ],start_year= 0, saveOutput=False):
        super(TBTimeseriesAnalyzer, self).__init__(site)
        self.filenames = [filename]
        self.channels = channels
        self.group_function = group_function
        self.filter_function = filter_function
        self.select_function = select_function
        self.plot_function   = plot_function
        self.saveOutput = saveOutput
        self.start_year =  start_year

    def filter(self, sim_metadata):
        return self.filter_function(sim_metadata)

    def validate_channels(self, keys):
        self.channels = [c for c in self.channels if c in keys] if self.channels else keys

    def get_channel_data(self, data_by_channel, header=None):
        channel_series = [self.select_function(data_by_channel[channel]) for channel in self.channels]
        return pd.concat(channel_series, axis=1, keys=self.channels)

    def apply(self, parser):
        data = parser.raw_data[self.filenames[0]]
        data_by_channel = data
        self.validate_channels(list(data.columns.values))
        channel_data = self.get_channel_data(data_by_channel, list([data.columns.values]))
        channel_data.group = self.group_function(parser.sim_id, parser.sim_data)
        channel_data.sim_id = parser.sim_id
        return channel_data

    def combine(self, parsers):
        # Gathering selected data from parser threads...
        selected = [p.selected_data[id(self)] for p in parsers.values() if id(self) in p.selected_data]

        # Combining selected data...
        combined = pd.concat(selected, axis=1, 
                             keys=[(d.group, d.sim_id) for d in selected], 
                             names=self.data_group_names)

        # Re-ordering multi-index levels...
        self.data = combined.reorder_levels(self.ordered_levels, axis=1).sort_index(axis=1)
        #self.data['channel']['Year']  = self.data['channel']['Year'] + self.start_year
        self.data['Year'] = self.data['Year'].apply(lambda y: y + self.start_year, axis=1)
        #self.data.set_index('Year', drop= False,inplace= True)

    def finalize(self):
        if self.saveOutput:
            self.data.to_csv(self.output_file)

    def plot(self):
        index_input =  self.data['Year']
        self.plot_channel_on_axes = lambda channel, ax: self.plot_function(self.data[channel].dropna() , ax, index_input)
        plot_by_channel(self.plot_name, self.channels, self.plot_channel_on_axes)
        import matplotlib.pyplot as plt
        plt.show()

