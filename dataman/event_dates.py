"""Module containing object for storing earthquake event dates, including 
uncertainties.

Jonathan Griffin
January 2019
"""

import numpy as np

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
        

        
