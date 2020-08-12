"""Module containing object for storing earthquake event dates, including 
uncertainties.

Jonathan Griffin
January 2019
"""

import numpy as np
import matplotlib
from matplotlib import pyplot
from matplotlib.patches import Ellipse
from scipy.stats import kde, spearmanr, gamma
from QuakeRates.utilities.inverse_transform_sample import ivt_expon
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
        if len(self.dates) > 1:
            bin_width = (self.dates[1] - self.dates[0])/2
            self.date_bins = list(self.dates - bin_width)
            self.date_bins.append(self.dates[-1] + bin_width)
        else:
            self.data_bins = self.dates[0]
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

    def __init__(self, event_list, name=None):
        """Intialise class
        :param event_list: List of EventDate objects ordered in forward
        running chronological order.
        """
        self.event_list = event_list
        self.num_events = len(event_list)
        self.name = name

    def gen_chronologies(self, n, search_limit=500, min_separation=20,
                         observation_start=None, observation_end=2019):
        """Generate n randomly sampled chronolgies for the events in
        EventSet object. As dating uncertainties may overlap bewteen events,
        event chronology is enforced and random samples that aren't in
        chronological order are discarded.
        :param n: Integer number of random samples to draw.
        :param min_separation: Integer number of minimum number of years that
        should separate consecutive events. It is assumed that some minimum
        period of time should have elasped for discrete palaeo-earthquakes
        to be observed.
        :param observation_start: Start of the observation period. If None,
        will be calculated as the data of the first earthquake in a sequence.
        But note that in some situations we may have other geological evidence
        of the start of the observation period.
        :param observation_end: Usually the present year, used for calculating
        the total length of the observation period.
        Ref: Biasi et al. 2002. Paleoseismic Event Dating and the Conditional
        Probability of Large Earthquakes on the Southern San Andreas Fault,
        California. Bulletin of the Seismological Society of America 92(7).
        """
        
        n_samples = 0
        n_tries = 0
        chron_samples = []
        loop_counter = 0
        self.add_events = False
        while n_samples < n:
            for i, event in enumerate(self.event_list):
                event.random_sample(n, plot=False)
                try:  # Append to previous sample if already exists
                    chron_samples[i] = chron_samples[i] + \
                        (event.random_from_cdf.tolist())
                except IndexError:
                    chron_samples.append(event.random_from_cdf.tolist())
            if loop_counter == 0:
                chron_tmp = np.vstack(chron_samples).T
                    #            chron_tmp = np.array(chron_samples).T
                # Now we add a future event as a random variable  using the conditional probability
                # assuming exponentially distributed inter-event times based
                # on mean of existing inter-event times
                # Only do if open interval is large (> mean ie time + 2 sigma).
                interevent_times = np.diff(chron_tmp.T, axis=0)
                #            print('interevent_times', interevent_times)
                mean_ie_time = np.mean(interevent_times)
                std_ie_time = np.std(interevent_times)
                #            print('mean_ie_time', mean_ie_time)
                #            print('chronologies.T[-1]', chronologies.T[-1])
                time_elapsed = observation_end - np.mean(chron_tmp.T[-1])
                #            print('mean last event')
                if (time_elapsed) > \
                   (mean_ie_time + 2*std_ie_time):
                    self.add_events = True
                    if self.name.startswith('Alpine'):
                        # Ignore records not complete until present
                        # i.e. Berryman 2012
                        self.add_events = False
                    if self.name.endswith('noadd'): # For testing not doing this
                        self.add_events = False
                    if self.name.endswith('nexytear'):
                        self.add_events = False
            if self.add_events:
                future_events = ivt_expon(1/mean_ie_time, a=time_elapsed,
                                          b=np.inf, n_samples=n)
                future_events += observation_end
                future_events -= time_elapsed # Previously this was double counted
                future_events = list(future_events)
                if loop_counter == 0:
                    chron_samples.append(future_events)
                else:
                    chron_samples[-1] = chron_samples[-1] + future_events
            # Now we check for chronological order
            chronologies = np.array(chron_samples).T
            c = chronologies[~np.any(np.diff(chronologies)<min_separation,
                                     axis=1)]
            chron_samples = c.T.tolist()
            n_samples = len(chron_samples[0])
            loop_counter += 1
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
        if observation_start is None:
            self.observation_period = np.abs(observation_end -  c[:,0])
        else:
            self.observation_period = np.abs(
                np.ones(n)* (observation_end - observation_start))
#        print('Length of observation periods is ',
#              self.observation_period, ' years.')
        self.chronology = c.T

    def plot_chronology(self, fig_filename, normalise=False):
        if hasattr(self, 'chronology'):
            pyplot.clf()
            if normalise:
                self.chronology_normalised = []
            for i, event in enumerate(self.event_list):
                if normalise:
                    # Normalise distirbutions by dividing by maximum value
                    event.probabilities_normalised = event.probabilities / \
                        max(event.probabilities)
                    self.chronology_normalised.append(self.chronology[i] / \
                                                      max(self.chronology[i]))
                    pyplot.plot(event.dates, event.probabilities_normalised,
                                color='k')
                    pyplot.hist(self.chronology_normalised[i],
                                bins=event.date_bins, density=True)
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
        self.interevent_times = np.diff(self.chronology, axis=0)
        self.means = np.mean(self.interevent_times, axis=0)
        self.stds = np.std(self.interevent_times, axis=0)
        print('Mean recurrence interval', np.mean(self.means))
        print('Recurrence interval standard devation',
              np.mean(self.stds))
        self.covs = self.stds/self.means
        # Now calculate normalised COV
        # (Brustiness parameter - Goh and Barabasi 2008; Chen et al 2020)
        self.burstiness = (self.stds - self.means)/ \
            (self.stds + self.means)
        self.mean_burstiness = np.mean(self.burstiness)
        self.std_burstiness = np.std(self.burstiness)
        self.mean_cov = np.mean(self.covs)
        print('Mean COV', self.mean_cov)
        print('Min COV', np.min(self.covs))
        print('Max COV', np.max(self.covs))
        # Now calculate bounds on COV and burstiness distribution
        self.cov_lb = np.percentile(self.covs, 2.5)
        self.cov_ub = np.percentile(self.covs, 97.5)
        self.burstiness_lb = np.percentile(self.burstiness, 2.5)
        self.burstiness_ub = np.percentile(self.burstiness, 97.5)
        pyplot.clf()
        pyplot.hist(self.covs, bins=25, density=True, edgecolor='0.2',
                    facecolor='0.6')
        pyplot.xlabel('Coefficient of variation')
        pyplot.ylabel('Probability density')
        pyplot.savefig('covs.png')
        # Plot burstiness distribution
        pyplot.clf()
        pyplot.hist(self.burstiness, bins=25, density=True, edgecolor='0.2',
                    facecolor='0.6')
        pyplot.xlabel('Burstiness')
        pyplot.ylabel('Probability density')
        pyplot.savefig('burstiness.png')

    def cov_density(self, fig_filename=None):
        """Calculate the density of the earthquake record's COV and
        long-term rate from the samples.
        """
        try:                                                                                                                                                                                         
            self.interevent_times
        except AttributeError as err:
            e = 'Need to call self.gen_chronologies before COV calculations'                                                           
            print(e)                                                                                                                  
            raise 
        self.long_term_rates = self.num_events/self.observation_period
                                                    
        """
        # Calculate error ellipses
        covariance = np.cov(self.covs, self.long_term_rates)
        vals, vecs = np.linalg.eig(covariance)
        order = vals.argsort()[::-1]
        vals = vals[order]
        vals = np.sqrt(vals)
        vecs = vecs[order]
        theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
        pyplot.clf()
        ax = pyplot.subplot(111)#, aspect='equal')
        ax.scatter(self.covs, self.long_term_rates, c='0.6', s=1)
        for j in range(1, 4):
            ell = Ellipse(xy=(np.mean(self.covs), np.mean(self.long_term_rates)),
                width=vals[0]*j*2, height=vals[1]*j*2,
                angle=theta)#np.rad2deg(np.arccos(v[0, 0])))
            ell.set_facecolor('none')
            ell.set_edgecolor('0.2')
            ax.add_artist(ell)
        """
        if fig_filename is not None:
            # Plot as contoured density plot
            nbins = 100
            cov_samples = np.array([self.covs, self.long_term_rates])
            x, y = cov_samples
            # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
            k = kde.gaussian_kde(cov_samples)
            # Slightly extend bounds for smoother plotting
            xi, yi = np.mgrid[0.98*x.min():1.02*x.max():nbins*1j, 0.98*y.min():1.02*y.max():nbins*1j]
            zi = k(np.vstack([xi.flatten(), yi.flatten()]))
            # Calculate percentiles
            cumsum = np.cumsum(zi)
            # Normalise data and calculate bottom 5th percentile for drawing contour
            # that contains 95% of the distribution.
            zi_norm = zi/max(cumsum)
            perc_5th = 0.05*max(zi_norm)
            # Plot the data
            pyplot.clf()
            ax = pyplot.subplot(111)
            # Slightly extend the bounds of the data for smoother plotting
            ax.pcolormesh(xi, yi, zi_norm.reshape(xi.shape), shading='gouraud', cmap=pyplot.cm.BuGn)
            ax.contour(xi, yi, zi_norm.reshape(xi.shape), [perc_5th])
            #FIXME - make bounds parametric/automatic
            ax.set_xlim([0,3])
            ax.set_ylim([1./100000, 1./100])
            ax.set_yscale('log')
            ax.set_xlabel('COV')
            ax.set_ylabel('Long-term rate (events per year)')
            pyplot.savefig(fig_filename)

    def basic_chronology_stats(self):
        """ Generate some basic statistics from the set of chronolgies
        """

        # Take global mean and standard devation of interevent times
        # across all simulated chronologies
        self.mean_interevent_time = np.mean(self.interevent_times)
        self.std_interevent_time = np.std(self.interevent_times)

        # Calculate mean value of shortest two consecutive interevent times.
        # This aims to estimate a 'short term rate'
        # Calculate sum of consecutive interevent_times
        min_interevent_times = []
        min_interevent_pairs = []
        max_interevent_times = []
        for int_ev_set in self.interevent_times.T:
            int_ev_pairs = np.array([np.sum(int_ev_set[current: current+2]) for current \
                                     in range(0, len(int_ev_set))])[:-1]
            min_int_ev_pair = np.min(int_ev_pairs)
            min_interevent_pairs.append(min_int_ev_pair)
            max_interevent_times.append(np.max(int_ev_set))
            min_interevent_times.append(np.min(int_ev_set)) 
        # Find mean and std of shortest two consecutive interevent times
        #self.minimum_two_interevent_times = np.mean(min_interevent_pairs)
        # Divide by two to estimate mean value of within cluster rate
        min_interevent_pairs = np.array(min_interevent_pairs)
        max_interevent_times = np.array(max_interevent_times)
        min_interevent_times = np.array(min_interevent_times)
        self.mean_minimum_pair_interevent_time =  np.mean(min_interevent_pairs/2.)
        self.std_minimum_pair_interevent_time =  np.std(min_interevent_pairs/2.)
        self.minimum_pair_interevent_time_lb = np.percentile((min_interevent_pairs/2.), 2.5)
        self.minimum_pair_interevent_time_ub = np.percentile((min_interevent_pairs/2.), 97.5)
        print('Minimum_pair_interevent time', self.mean_minimum_pair_interevent_time)
        # Get mean and std of maximum interevent times
        self.mean_maximum_interevent_time = np.mean(max_interevent_times)
        self.std_maximum_interevent_time = np.std(max_interevent_times)
        self.maximum_interevent_time_lb = np.percentile((max_interevent_times/2.), 2.5)
        self.maximum_interevent_time_ub = np.percentile((max_interevent_times/2.), 97.5)
        print('self.mean_maximum_interevent_time', self.mean_maximum_interevent_time)
        # Get mean and std of minimum interevent times
        self.mean_minimum_interevent_time = np.mean(min_interevent_times)
        self.std_minimum_interevent_time = np.std(min_interevent_times)
        self.minimum_interevent_time_lb = np.percentile((min_interevent_times/2.), 2.5)
        self.minimum_interevent_time_ub = np.percentile((min_interevent_times/2.), 97.5)
        print('self.mean_minimum_interevent_time', self.mean_minimum_interevent_time)
        # Get ratios
        self.ratio_min_pair_max = min_interevent_pairs/max_interevent_times
        self.ratio_min_max = min_interevent_times/max_interevent_times
        self.mean_ratio_min_pair_max = np.mean(self.ratio_min_pair_max)
        self.mean_ratio_min_max = np.mean(self.ratio_min_max)
        self.std_ratio_min_pair_max = np.std(self.ratio_min_pair_max)
        self.std_ratio_min_max = np.std(self.ratio_min_max)
        self.ratio_min_pair_max_lb = np.percentile(self.ratio_min_pair_max, 2.5)
        self.ratio_min_pair_max_ub = np.percentile(self.ratio_min_pair_max, 97.5)
        self.ratio_min_max_lb = np.percentile(self.ratio_min_max, 2.5)
        self.ratio_min_max_ub = np.percentile(self.ratio_min_max, 97.5) 

    def memory_coefficient(self):
        """ Calculate the memory coeeficient of 
        Goh, K.-I. and A.-L. Barabási (2008). Burstiness and memory in 
        complex systems, Europhysics Lett., 81, no. 4, 48002.
        """
        #print(self.interevent_times)
        self.m1 = np.mean(self.interevent_times[0:-1], axis=0)
        self.m2 = np.mean(self.interevent_times[1:], axis = 0)
        self.s1 = np.std(self.interevent_times[0:-1], axis=0)
        self.s2 = np.std(self.interevent_times[1:], axis = 0)
        # Handle case for small number of interevent times (e.g. 2) where we
        # can randomly generate identifical values and get a std=0
        # Give this a small non-zero value to avoid divide by zero
        self.s1[np.argwhere(self.s1==0)] = 1
        self.s2[np.argwhere(self.s2==0)] = 1
        self.s1s2 = self.s1*self.s2
        # Loop over chronologies
        ie_min_m1 = self.interevent_times - self.m1
        # Remove last element
        ie_min_m1 = ie_min_m1[0:-1]
        ie_min_m2 = self.interevent_times - self.m2
        # Remove first element
        ie_min_m2 = ie_min_m2[1:]
        numerator = ie_min_m1 * ie_min_m2
        numerator = np.sum(numerator, axis=0)
        sum_term = numerator/self.s1s2
        self.mem_coef = sum_term * (1/(self.num_events - 1))
        self.mean_mem_coef = np.mean(self.mem_coef)
        print('Mean memory coefficient', self.mean_mem_coef)
        self.memory_lb = np.percentile(self.mem_coef, 2.5)
        self.memory_ub = np.percentile(self.mem_coef, 97.5)

    def memory_spearman_rank_correlation(self):
        """ Calculate alternative memory coefficient using Spearman rank
        correlation. This may avoid some biases in the calculation. 
        See: Schleiss, M. and J. A. Smith (2016). Two Simple Metrics for 
        Quantifying Rainfall Intermittency: The Burstiness and Memory of
        Interamount Times, J. Hydrometeorol., 17, no. 1, 421–436.
        """
        # Lag-1
        a = self.interevent_times[0:-1].T
        b = self.interevent_times[1:].T
        # Lag-2
        a2 = self.interevent_times[0:-2].T
        b2 = self.interevent_times[2:].T 
#        print(a, b)
        rhos = []
        rhos2 = []
        #        self.rhos, self.pvalues = spearmanr(a, b)
        for i, ie_t in enumerate(a):
#            print(ie_t.T)
#            print(b[i].T)
            rho, pvalue = spearmanr(ie_t.T, b[i].T)
#            print('rho', rho)
            rhos.append(rho)
        for i, ie_t in enumerate(a2):
            rho, pvalue = spearmanr(ie_t.T, b2[i].T)
            rhos2.append(rho)
        self.rhos = np.array(rhos)
        self.rhos2 = np.array(rhos2)
#        print(self.rhos, self.pvalues)
        self.mean_rho = np.mean(self.rhos)
        self.mean_rho2 = np.mean(self.rhos2)
        print('Mean Spearman Rank coefficient', self.mean_rho)
        print('Mean Spearman Rank coefficient', self.mean_rho2)
        self.rho_lb = np.percentile(self.rhos, 2.5)
        self.rho_ub = np.percentile(self.rhos, 97.5)
        self.rho2_lb = np.percentile(self.rhos2, 2.5)
        self.rho2_ub = np.percentile(self.rhos2, 97.5) 
#        print(self.rho_lb, self.rho_ub)
    
    def plot_interevent_time_hist(self, fig_filename='interevent_times.png'):
        """ Plot histogram of interevent times
        """
        ie_times_flat = self.interevent_times.flatten()
        max_interevent_time = max(ie_times_flat)
        if max_interevent_time < 200:
            s = 10
        elif max_interevent_time < 1000:
            s = 50
        elif max_interevent_time < 10000:
            s = 200
        elif max_interevent_time < 1e5:
            s = 2000
        else:
            s = 5000
        bins = np.arange(1, max_interevent_time, s)
        pyplot.clf()
        pyplot.hist(ie_times_flat, bins = 50,
                    facecolor='0.6', edgecolor='0.2', density=True)
        pyplot.savefig(fig_filename)

    def fit_gamma(self):
        """ Fit a gamma distribution to each chronology, get mean stats
        """
        self.gamma_alphas = []
        for ie_set in self.interevent_times.T:
            gamfit = gamma.fit(ie_set, floc=0)
            self.gamma_alphas.append(gamfit[0])
        # test fitting to all data at once
        gamfit_all = gamma.fit(self.interevent_times.flatten(), floc=0)
        self.gamma_alphas = np.array(self.gamma_alphas)
        self.mean_gamma_alpha = np.mean(self.gamma_alphas)
        self.mean_gamma_alpha_all = gamfit_all[0]        
