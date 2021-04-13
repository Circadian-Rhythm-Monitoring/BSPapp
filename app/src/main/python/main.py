import mat4py as loadmat
import numpy as np
from numpy import random

import scipy as sci
from scipy import signal
from scipy.fft import fft, ifft
from scipy.special import comb

import math as math
import scipy.io

from GainOpt_FilterDyn_Class import GainOpt_FilterDyn


def main():

    pi = math.pi
    omg = (2 * pi) / 24
    zeta = 1
    gamma_d = 1
    gamma_omg = 0
    order = 1
    stateLength = (2 * order + 1)
    subject = 3

    newFilter = GainOpt_FilterDyn(subject,order,stateLength,omg,zeta,gamma_d,gamma_omg)
    final = newFilter.Final
    return(final)