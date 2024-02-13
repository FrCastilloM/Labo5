    # -*- coding: utf-8 -*-
"""
Generador de funciones Tektronix AFG 3021B
Manual U (web): https://github.com/diegoshalom/labosdf/blob/master/manuales/AFG3021B%20user%20manual.pdf
Manual P (web): https://github.com/diegoshalom/labosdf/blob/master/manuales/AFG3021B%20Programmer%20Manual.pdf
"""


import time

import numpy as np
import pyvisa as visa

class AFG3021B:
    
    def __init__(self, name='USB0::0x0699::0x0346::C034165::INSTR'):
        self._generador = visa.ResourceManager().open_resource(name)
        print(self._generador.query('*IDN?'))
        
        #Activa la salida
        self._generador.write('OUTPut1:STATe on')
        # self.setFrequency(1000)
        
    # def __del__(self):
    #     self._generador.close()
        
    def setFrequency(self, freq):
        self._generador.write(f'FREQ {freq}')
        
    def getFrequency(self):
        return self._generador.query_ascii_values('FREQ?')
    
    def setOffset(self, offset, channel = 1): #gen.SetOffset(1) V
        self._generador.write('SOURce{}:VOLTage:LEVel:IMMediate:OFFSet {}'.format(channel, offset))
        
    def setAmplitude(self, freq):
        print('falta')
        
    def getAmplitude(self):
        print('falta')
        return 0



# #generador de funciones
# fungen = AFG3021B(name = 'USB0::0x0699::0x0346::C034198::INSTR')
# fungen.getFrequency()

# #barrido de frecuencia
# for freq in range(1000,5000,1000):
#     print(freq)
#     fungen.setFrequency(freq)
#     time.sleep(0.1)
