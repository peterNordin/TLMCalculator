import numpy as np

__author__ = "Peter Nordin"
__license__ = "GPLv3"
__email__ = "peter.nordin@liu.se"
'''This module contains basic TLM classes for model approximation calculations'''


def notzero(value, eps):
    if abs(value) < eps:
        if value < 0:
            return -eps
        else:
            return eps
    return value


class TLMParams:
    def __init__(self):
        self.T = 0
        self.zc = 0
        self.length = 0
        self.wavespeed = 0
        self.ind = 0
        self.stiffness = 0
        self.cap = 0
        self.ind2 = 0
        self.cap2 = 0


class TLMParasiticParams:
    def __init__(self):
        self.ind = 0
        self.cap = 0


class TLMLineParameters:
    def __init__(self, wavespeed, length, zc, mass, capacitance, zctype='normal' ):
        self.desired = TLMParams()
        self.desired.wavespeed = wavespeed
        self.desired.T = length/wavespeed
        self.desired.length = length
        self.desired.zc = zc
        self.desired.ind = mass
        self.desired.stiffness = 1/capacitance
        self.desired.cap = capacitance

        self.parasitic = TLMParasiticParams()
        self.parasitic.ind = 0
        self.parasitic.cap = 0

        self.percivedLC = TLMParasiticParams()

        self.percived = TLMParams()
        self.percived.T = 0
        self.percived.zc = 0
        self.percived.length = 0
        self.percived.ind = 0
        self.percived.stiffness = 0
        self.percived.cap = 0

        self.zctype = zctype

    def calcStuff(self, ts, n, overrideZc=None):
        # print('in tlm calc stuff')

        # Set n automatically
        if n == -1:
            n = max(round(self.desired.T/ts), 1)
            print('n: '+str(n))

        T = n*ts                                           # Perceived TLM delay, numStep*sampletime

        if overrideZc is None:
            zc = self.calcZc(T)
        else:
            zc = overrideZc

        self.percived.zc = zc
        self.percived.T = T
        self.percived.length = T*self.desired.wavespeed    # Perceived dealy * wavespeed

        self.percived.ind = T*zc
        self.percived.cap = T/zc
        self.parasitic.ind = self.percived.ind - self.desired.ind
        self.parasitic.cap = self.percived.cap - self.desired.cap

        if self.zctype == 'normal':
            self.percivedLC.ind = T*np.sqrt(self.desired.ind/self.desired.cap)
            self.percivedLC.cap = T/np.sqrt(self.desired.ind/self.desired.cap)
        elif self.zctype == 'pureinductance':
            self.percivedLC.ind = T*self.desired.ind/T
            self.percivedLC.cap = T*T/self.desired.ind
        elif self.zctype == 'purecapacitance':
            self.percivedLC.ind = T*self.desired.cap/T
            self.percivedLC.cap = T*T/self.desired.cap
        else:
            raise Exception('Wrong zctype')

        self.percived.stiffness = 1/notzero(self.percived.cap, 1e-100)

    def calcZc(self, T):
        if self.zctype == 'normal':
            zc = np.sqrt(self.desired.ind/self.desired.cap)
        elif self.zctype == 'pureinductance':
            zc = self.desired.ind/T
        elif self.zctype == 'purecapacitance':
            zc = self.desired.stiffness*T
        else:
            raise Exception('Wrong zctype')
        # print('zc: '+str(zc))
        return zc


class TLMHydraulicLine(TLMLineParameters):
    def __init__(self, bulk, rho, length, area, zctype):
        self.volume = area*length
        self.bulk = bulk
        waveSpeed = np.sqrt(bulk/rho)
        ind = self.volume*rho
        cap = self.volume/self.bulk

        TLMLineParameters.__init__(self, waveSpeed, length, np.sqrt(ind/cap), ind, cap)
        self.zctype = zctype


class TLMMechanicLine(TLMLineParameters):
    def __init__(self, young, rho, length, area, zctype):
        self.volume = area*length
        self.youngsmodulus = young
        wavespeed = np.sqrt(young/rho)

        # print(area)
        desired_inductance = self.volume*rho
        desired_stiffness = (young*area)/length
        desired_zc = np.sqrt(desired_inductance*desired_stiffness)

        TLMLineParameters.__init__(self, wavespeed, length, desired_zc, desired_inductance, 1/desired_stiffness)
        self.zctype = zctype


class TLMMechanicSpring(TLMLineParameters):
    def __init__(self, desired_stiffness, T):
        desired_inductance = 0
        desired_zc = desired_stiffness * T

        TLMLineParameters.__init__(self, 1e99, 0, desired_zc, desired_inductance, 1/desired_stiffness)
        self.zctype = 'purecapacitance'


class TLMHydraulicSpring(TLMLineParameters):
    def __init__(self, bulk, volume, T):
        desired_inductance = 0
        desired_zc = bulk/volume * T
        desired_cap = volume/bulk

        TLMLineParameters.__init__(self, 1e99, 0, desired_zc, desired_inductance, desired_cap)
        self.zctype = 'purecapacitance'
