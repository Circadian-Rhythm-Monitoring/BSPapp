import mat4py as loadmat
import numpy as np
from numpy import random
np.random.seed(1)

import scipy as sci
from scipy import signal
from scipy.fft import fft, ifft
from scipy.special import comb

from os.path import dirname, join

import math as math
import scipy.io
import itertools as itertools


class GainOpt_FilterDyn:
    # Optimization settings
    __max_iterations = 10  # 15  #Max iterations
    __rho = 2

    #GLOBAL VARIABLESsss
    __mu = 50 #100
    __lambda_ = 10 #50
    Cost = np.zeros(__mu)

    #initlial conditions
    __pi = math.pi
    __omg = (2 * __pi) / 24
    __zeta = 1
    __gamma_d = 1
    __gamma_omg = 0
    __order = 1
    __stateLength = (2 * __order + 1)
    __subject = 3
    __len_ = 0  # Length of signal

    # ----------------Generating Initial __population
    __LB = -3  # 1*(10**-20)
    __UB = 20  # .1*(10**-2)

    __lStart = 0
    __rEnd = 4

    __len_ = 0

    #----------------Publicaly acsessalbe elements
    avgCost = None
    Ad = None
    Bd = None
    Cd = None
    Dd = None
    Final = None
    __combinations = None

    def __init__ (self,subject_,order_,stateLength_,omg_,zeta_,gamma_d,gamma_omg):

        __omg = omg_# (2 * __pi) / 24
        __zeta =  zeta_#1
        __gamma_d = gamma_d#1
        __gamma_omg = gamma_omg#0
        __order = order_#1
        __stateLength = stateLength_#(2 * __order + 1)
        __subject = subject_#3

        mid = (self.__rEnd + self.__lStart) / 2
        N = self.__rEnd * np.random.rand(self.__mu, self.__stateLength)
        zeroSize = (len(N), len(N[0]))
        self.__population = np.zeros(zeroSize)
        for i in range(0, len(N)):
            for j in range(0, len(N[0])):
                if N[i][j] >= mid:
                    self.__population[i][j] = 10 ** ((N[i][j]) - mid + self.__LB)
                else:
                    self.__population[i][j] = -10 ** (mid - (N[i][j]) + self.__LB)

        #-----------------Load the data in
        self.__t, self.__y = self.__loadData(self.__subject)

        # ----------------Fourier Transform of original signal
        # Take the Fourier Transform of the original signal to use in the
        # calculation of the costs
        T = self.__t[2] - self.__t[1]  # Sampling period
        Fs = 1 / T  # Sampling frequency
        self.__len_ = len(self.__t)  # Length of signal
        Y = fft(self.__y)  # True output to compare with ANF output
        P2 = abs(Y / self.__len_)  # Compute the two-sided spectrum P2 (Folded across the __y-axis)
        P = P2[0:math.floor(self.__len_ / 2)]  # Compute the single-sided spectrum P
        P[1:len(P) - 1] = 2 * P[1:len(P) - 1]  # Double everything except the DC term
        # f = np.zeros(math.floor(__len_/2))
        # for i in range(0,math.floor(__len_/2)):
        #     f[i] = Fs*(i)/__len_ # Define the freq based on the sampling rate
        f = np.array([Fs * i / self.__len_ for i in range(math.floor(self.__len_ / 2))])

        self.__combinations = list(itertools.combinations(range(1, self.__lambda_ + 1), 2))
        self.__combinations = np.matrix(self.__combinations)

        # should be creating a array of size __combinations x 1
        Labels = np.random.randint(1, len(self.__combinations), len(self.__combinations))

        # __generationCalc(Ad,Bd,Cd,__t,__y,P,f,Labels,False)
        self.__generationCalc(self.__t, self.__y, P, f, Labels, False)

        # This code is responsible for the gene pool filtering
        #   -sorts the cost array from lowest to highest
        #   -from the cost array, delete the top 50, do same for __population
        #   -feeds back into the __generationCalc function
        #   -__population is a global varaible and is updated as such
        self.avgCost = np.zeros(self.__max_iterations, dtype=float)
        for i in range(0, self.__max_iterations):
            self.__generationCalc(self.__t, self.__y, P, f, Labels, True)

            # This grabs the maxmum values of the childen and returns there index
            maxIndex = np.argpartition(self.Cost, -self.__lambda_)[-self.__lambda_:]

            # will remove the elements attached to the index
            self.Cost = np.delete(self.Cost, maxIndex)
            self.__population = np.delete(self.__population, maxIndex, axis=0)
            self.avgCost[i] = np.mean(self.Cost)

        self.Final = self.avgCost[len(self.avgCost) - 1]

    def __computeCost(self, yHat, P, f, len_):

        #to make sure yHat is a numpy array
        yHat = np.array(yHat)
        f = np.array(f)
        P1 = np.ones(math.ceil(len_/2))
        #computes the fft
        #   -Takes half of the signal and essentially flips it on itself
        #       This amounts to essentially just halfing the P2 array and then doubling it
        #       We can do this because the signal should be symetric
        Y = fft(yHat)
        Y = np.array(Y)
        P2 = np.absolute(Y / self.__len_)
        P1 = P2[0,0:math.ceil((self.__len_ / 2))]
        P1[1:len(P1)] = 2 * P1[1:len(P1)]
        P1= np.array(P1).T


        #This is the filtering part
        #   -Essentially we are cutting out frequencies we don'__t want
        #   -First find absolute value, filtering the unwanted frequencies
        #   -Second finds index where the value occured, then converts tuple to an int
        f_N1 = np.absolute(f-(1/24))
        f_N2 = np.absolute(f-.0289)
        N1 = np.where(f_N1 ==np.amin(f_N1))[0]
        N2 = np.where(f_N2 ==np.amin(f_N2))[0]
        NN = N1[0] - N2[0]

        P1_a = np.array(P1[0:NN])
        P1_b = np.array(P1[N1[0]-NN:N1[0]+NN+1])
        P_a = np.array(P[:NN])
        P_b = np.array(P[N1[0]-NN:N1[0]+NN+1])

        #Calculates the square error within the band each specified harmonic
        #harmonic and the DC term
        J_harm = np.trapz(np.square(P1_a - P_a)) + np.trapz(np.square(P1_b - P_b))

        P1_c = np.square(P1[NN :N1[0] - NN])
        P1_d = np.square(P1[N1[0] + NN + 1: len(P1)])

        #square of the signal outside the bands around each harmonic and beyond the last one
        J_noise = np.trapz(P1_c) + np.trapz(P1_d)

        # JHN = [J_harm,J_noise]
        return J_harm, J_noise

    def simulateDynamics(self, stateLength, t, y, Ad, Bd, Cd):
        xHat = np.zeros([stateLength,len(t)], dtype = float)
        xHat[:,0] = 1

        # first loop through columns,start one column in
        for j in range(1,len(t)):#columns
            #CK CODE FOR SIMULATE DYNAMICS

            xHat[:,j] = np.reshape(np.reshape(np.matmul(Ad,xHat[:,j-1]).T,(stateLength,1)) + (Bd*y[j-1]), (stateLength))

        yHat = np.matmul(Cd, xHat)

        #this is a matrix
        return [xHat,yHat]

    def __loadData(self, subject):
        self.__t = []
        self.__y = []
        tPlusY = []
        subjectData = []

        filename = join(dirname(__file__), "A3.mat")
        data = scipy.io.loadmat(filename)
        # data = scipy.io.loadmat('A' + str(subject) + '.mat')
        # takes the data from dictionary "data"
        subjectData = data['Data']

        i = 0
        # loop goes until the last row of the matrix
        while i < len(subjectData):
            self.__t.append(subjectData[i][0])
            self.__y.append(subjectData[i][len(subjectData[0]) - 1])
            i = i + 1
        self.__t = np.transpose(self.__t)

        # converting to numpy arrays for performance
        self.__t = np.array(self.__t)
        self.__y = np.array(self.__y)

        tPlusY = [self.__t, self.__y]

        #
        # data = scipy.io.loadmat('A' + str(__subject)+ '.mat')
        # #takes the data from dictionary "data"
        # subjectData = data['Data']
        # __t = subjectData[:,0]
        # __y = subjectData[:,1]

        return self.__t, self.__y


    def __createStateSpace (self, order, stateLength, t, omg, zeta, gamma_d, gamma_omg, L):
        Ac = np.zeros([stateLength,stateLength], dtype = float)
        Bc = np.zeros([stateLength,1],dtype = float)
        Cc = np.zeros([1,stateLength],dtype = float)
        Dc = 0

        #This is to populate State matrices
        #   -Python takes the LHS which forms a 2x2 matrix, matches it to the RHS which, to python looks like a
        #    2x2 matrix and assigns each element accordingly. Process follows suite for Bc matrix
        for k in range(order):
            i = (1 + k) * 2
            k = k + 1 #handles the one off error
            Ac[i - 2:i, i - 2:i] = [[0, 1], [float(-(k * omg) ** 2), 0]]
            Bc[i-2:i] = [[0],[(k*omg)**2]]
            Cc[0][i-1] = (2 * zeta) / (k * omg)

        #Assigns edge values
        Bc[len(Bc)-1][0] = gamma_d
        Cc[0][len(Cc[0])-1] = 1

        A = (Ac - L*Cc)
        B = L
        C = Cc
        D = Dc

        #uses the "scipy" library and converts the continuous matrixes into a SS
        #   -Then takes the continuous system and decritizes it (NOTE: method for c2d may not match model)
        dt = float(t[1]-t[0])
        contSystem = signal.StateSpace(A, B, C, D)
        discSystem = contSystem.to_discrete(dt)

        #Posible solution if .to_discrete system does not work.
        #   -function for computing matrix exponenetal (Future implimentation)
        #discSystem = signal.cont2discrete(contSystem,dt,'impulse')

        return discSystem

    # this function is design to handle both parts of the gene pool selection proccess
    #   -By setting the last argument to true, the process will be run on the additional elements
    #   -By setting the last argument to false, the process will be run to initialize __population
    #   -Cost and __population are global variables

    # __generationCalc(mat Ad, mat Bd, mat Cd, array __t, array __y, array __population, array P,array f) - return cost, __population, xHat, yHat
    # def __generationCalc(Ad,Bd,Cd, __t, __y, P, f,Labels, fExtra):
    def __generationCalc(self, t, y, P, f, Labels, fExtra):
        global population
        global Cost
        if fExtra == True:
            rangeCalc = self.__lambda_
        else:
            rangeCalc = self.__mu
        for i in range(0, rangeCalc):

            # this computes the remaining rows of the __population matrix
            if fExtra == True:
                k = self.__mu + i
                # this has to be checked.
                # test = np.mean(__population[__combinations[Labels[i],:],:],axis = 1)
                self.__population = np.append(self.__population, np.mean(self.__population[self.__combinations[Labels[i], :], :], axis=1), axis=0)

            # NOT TESTED
            # L = np.matrix(__population[i,:]).T

            # FORCING L FOR TESTING
            L = np.matrix([[0], [.0086], [.0339]])

            # NEWLY ADDED CODE FROM CK

            discSystem = self.__createStateSpace(self.__order, self.__stateLength, self.__t, self.__omg, self.__zeta, self.__gamma_d, self.__gamma_omg, L)

            self.Ad = discSystem.A
            self.Bd = discSystem.B
            self.Cd = discSystem.C
            self.Dd = discSystem.D

            eigenvalues = np.linalg.eig(self.Ad)

            # Catches errors
            if np.sum(np.absolute(eigenvalues[1]) > 1) > 0:
                if fExtra == True:
                    self.Cost[k] = intmax
                else:
                    self.Cost[i] = intmax
                continue

            # runs the system dynamics as well as the cost function
            # xHatyHat = simulateDynamics(__stateLength,__t,__y,Ad,Bd,Cd,L)
            xHatyHat = self.simulateDynamics(self.__stateLength, self.__t, self.__y, self.Ad, self.Bd, self.Cd)
            xHat = xHatyHat[0]
            yHat = xHatyHat[1]
            J_harm, J_noise = self.__computeCost(yHat, P, f, self.__len_)
            if fExtra == True:
                self.Cost = np.append(self.Cost, J_harm + J_noise)
                # Cost[k] = JHJN[0]+JHJN[1]
            else:
                self.Cost[i] = J_harm + J_noise