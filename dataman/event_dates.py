"""Module containing object for storing earthquake event dates, including 
uncertainties.

Jonathan Griffin
January 2019
"""

import numpy as np
import matplotlib
from matplotlib import pyplot

matplotlib.use('Agg')
np.random.seed(23)


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
        pyplot.plot(self.dates, self.probabilities, color='k')
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
        self.random_from_cdf = self.dates[value_bins]
        bin_width = (self.dates[1] - self.dates[0])/2
        self.date_bins = list(self.dates - bin_width)
        self.date_bins.append(self.dates[-1] + bin_width)
        if plot:
            pyplot.clf()
            pyplot.plot(self.dates, self.probabilities, color='k')
            pyplot.hist(self.random_from_cdf, bins=self.date_bins,
                        facecolor='0.6', edgecolor='0.2', density=True)
            pyplot.savefig(fig_filename)
#        self.random_from_cdf

class EventSet(object):
    """Class for storing an ordered list of earthquake events.
    This allows random samples of the event chronology to be drawn
    """

    def __init__(self, event_list):
        """Intialise class
        :param event_list: List of EventDate objects ordered in forward
        running chronological order.
        """
        self.event_list = event_list

    def gen_chronologies(self, n, search_limit=10, min_separation=20):
        """Generate n randomly sampled chronolgies for the events in
        EventSet object. As dating uncertainties may overlap bewteen events,
        event chronology is enforced and random samples that aren't in
        chronological order are discarded.
        :param n: Integer number of random samples to draw.
        :param min_separation: Integer number of minimum number of years that
        should separate consecutive events. It is assumed that some minimum
        period of time should have elasped for discrete palaeo-earthquakes
        to be observed.
        Ref: Biasi et al. 2002. Paleoseismic Event Dating and the Conditional
        Probability of Large Earthquakes on the Southern San Andreas Fault,
        California. Bulletin of the Seismological Society of America 92(7).
        """
        
        n_samples = 0
        n_tries = 0
        chron_samples = []
        while n_samples < n:
            for i, event in enumerate(self.event_list):
                event.random_sample(n, plot=False)
                try:  # Append to previous sample if already exists
                    chron_samples[i] = chron_samples[i] + \
                        (event.random_from_cdf.tolist())
                except IndexError:
                    chron_samples.append(event.random_from_cdf.tolist())
            # Now we check for chronological order
            chronologies = np.array(chron_samples).T
            c = chronologies[~np.any(np.diff(chronologies)<min_separation,
                                     axis=1)]
            chron_samples = c.T.tolist()
            n_samples = len(chron_samples[0])
#            print(n_samples)
            n_tries += n

            msg = 'Could not find ' + str(n) + ' samples in ' + \
                'chronological order after drawing ' + \
                str(search_limit*n) + ' random samples.' + \
                'Please  check the input event order or increase ' + \
                'the search_limit parameter.'
            assert n_tries < search_limit*n, msg
        # Now need to clip to only have n samples, if more than n generated.
        c = c[0:(n)]
        print('Number of chronology samples', c[:,0].size)
        self.chronology = c.T

    def plot_chronology(self, fig_filename, normalise=False):
        if hasattr(self, 'chronology'):
            pyplot.clf()
            if normalise:
                self.chronology_normalised=[]
            for i, event in enumerate(self.event_list):
                if normalise:
                    # Normalise distirbutions by dividing by maximum value
                    event.probabilities_normalised = event.probabilities / \
                        max(event.probabilities)
                    self.chronology_normalised.append(self.chronology[i] / \
                                                      max(self.chronology[i]))
                    pyplot.plot(event.dates, event.probabilities_normalised, color='k')
                    pyplot.hist(self.chronology_normalised[i], bins=event.date_bins,                                            
                            density=True)

                else:
                    pyplot.plot(event.dates, event.probabilities, color='k')
                    pyplot.hist(self.chronology[i], bins=event.date_bins,
                                density=True)#, edgecolor='0.2')
                    # facecolor='0.6', edgecolor='0.2', density=True)
                ax = pyplot.gca()
                ax.set_xlabel('Years')
                ax.set_ylabel('Probability')
                #ax.set_yscale('log')
            pyplot.savefig(fig_filename)
        else:
            print('Need to call self.gen_chronologies before plot_chronology')

    def write_chronology(self, filename):
        """Dump chronology samples to .csv file
        """
        if hasattr(self, 'chronology'):
            np.savetxt(filename, self.chronology.T, delimiter=',')
        else:
            print('Need to call self.gen_chronologies before write_chronology')
            
        
    def calculate_cov(self):
        """Calculate the coeeficient of variation from the randomly 
        sampled chronolgies.
        """
        try:
            self.chronology
        except AttributeError as err:
            e = 'Need to call self.gen_chronologies before COV calculations' 
            print(e)
            raise 
        # The nth row of interevent_times contains all reaslisations of
        # the interevent time between the nth and n+1 event (counting
        # forward in time)
        interevent_times = np.diff(self.chronology, axis=0)
#        print(self.chronology, len(self.chronology))
#        print(interevent_times, len(interevent_times))
        means = np.mean(interevent_times, axis=0)
        stds = np.std(interevent_times, axis=0)
        print('Mean recurrence interval', np.mean(means))
        print('Recurrence interval standard devation',
              np.mean(stds))
#        print(stds)
        self.covs = stds/means
#        print(self.covs)
        print('Mean COV', np.mean(self.covs))
        print('Min COV', np.min(self.covs))
        print('Max COV', np.max(self.covs))
        pyplot.clf()
        pyplot.hist(self.covs, bins=25, density=True, edgecolor='0.2',
                    facecolor='0.6')
        pyplot.xlabel('Coefficient of variation')
        pyplot.ylabel('Probability density')
        pyplot.savefig('covs.png')
        
