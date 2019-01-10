"""Module containing object for storing earthquake event dates, including 
uncertainties.

Jonathan Griffin
January 2019
"""

import numpy as np
import matplotlib
from matplotlib import pyplot
matplotlib.use('Agg')

class EventDate(object):
    """Class for storing informationg related to a single earthquake event
    """

    def __init__(self, id, op, date_type):
        """
        :params id: ID for the earthquake event
        :params date_op: string, equivalent to OxCal 'op' parameter, i.e.
        whether the date is a raw date '(R_14) or a calculated event date
        ('calculate').
        :params date_type: string, equivalent to OxCal 'type' parameter.
        Can be either 'likelihood' or 'posterior'.
        """

        self.id = id
        self.date_type = date_type

    def add_dates_and_probs(self, dates, probabilities):
        """Adds list of dates and their associated probabilities
        to the event
        :params dates: Array or list (Nx1) of dates associated with the event
        :params probabilities: Array or list (Nx1) of probabilities associated
        with each of the dates in dates.
        """

        self.dates = np.array(dates)
        self.probabilities = np.array(probabilities)
        
    def plot_date_pdf(self, fig_filename='event.pdf'):
        """Produce a simple plot of the probability distribution for the event 
        date
        """

        pyplot.clf()
        pyplot.plot(self.dates, self.probabilities)
        pyplot.xlabel('Date')
        pyplot.ylabel('Probability')
        pyplot.title(self.id)
        pyplot.savefig(fig_filename)
        
    def random_sample(self, n, plot=False, fig_filename='random_sample.png'):
        """Generate n random samples from the probability distribution
        of the dates.
        :param n: integer, number of random samples to draw
        :param plot: If True, plot the random samples
        """
        cdf = np.cumsum(self.probabilities)
        cdf = cdf / cdf[-1]
        values = np.random.rand(n)
        value_bins = np.searchsorted(cdf, values)
        random_from_cdf = self.dates[value_bins]

        if plot:
            pyplot.clf()
            pyplot.plot(self.dates, self.probabilities)
            pyplot.hist(random_from_cdf, len(self.dates))
            #pyplot.plot(value_bins, random_from_cdf)
            pyplot.savefig(fig_filename)
