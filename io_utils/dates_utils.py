from scipy.ndimage import binary_fill_holes

import math
import numpy as np

from datetime import date, datetime, timedelta

def get_days_from_month(month, year=datetime.now().year):
    """
    Gets a list of integers with the days of the month and day of the year for and specific month and year
    :param month:
    :param year:
    :return:
    """
    first_jan = date(year, 1, 1)
    first_month = date(year, month, 1)
    if month == 12:
        next_month = date(year+1, 1, 1)
    else:
        next_month = date(year, month+1, 1)

    days_of_year = np.arange(first_month.toordinal()-first_jan.toordinal()+1, next_month.toordinal()-first_jan.toordinal()+1)
    days_of_month = np.arange(0,  next_month.toordinal() - first_month.toordinal())
    return days_of_month, days_of_year  # TODO compute from the month

def get_day_of_year_from_month_and_day(month, day_of_month, year=datetime.now().year):
    """
    Gets a list of integers with the days of the month, starting from 0 and from the day of the year
    :param month:
    :param year:
    :return:
    """
    first_jan = date(year, 1, 1)
    day_of_year = date(year, month, day_of_month).toordinal() - first_jan.toordinal() + 1
    return day_of_year

def get_month_and_day_of_month_from_day_of_year(day_of_year, year=datetime.now().year):
    """
    Gets a list of integers with the days of the month, starting from 0 and from the day of the year
    :param month:
    :param year:
    :return:
    """
    date_jan = date(year, 1, 1)
    curr_date = date_jan + timedelta(days=int(day_of_year)-1)
    return curr_date.month, curr_date.day
