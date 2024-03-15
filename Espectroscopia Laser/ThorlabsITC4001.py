# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 09:55:38 2024

@author: Franc
"""

import pyvisa as visa

class ITC4001:
    
    def __init__(self, name):
        self._thor = visa.ResourceManager().open_resource(name)
        print(self._thor.query("*IDN?"))
    	
    def __del__(self):
        self._thor.close()			

    def setCurrent(self, current):
        self._thor.write('SOUR:CURR {}'.format(current))
    
    def getCurrent(self):
        return self._thor.query('MEAS:CURR?')
    
    # def setTemp(self, temp):
        # self._thor.write('SOUR:TEMP[:SPO] {}C'.format(temp)) # NO FUNCIONA
    
    def getTemp(self):
        return self._thor.query('MEAS:TEMP?')
    
    def setModulation(self, freq, mode = 'TRI'): # mode = 'SIN', 'SQU', 'TRI'
        self._thor.write('SOUR:AM[:STAT] ON')
        self._thor.write('SOUR:AM:INT:SHAP {}'.format(mode))
        self._thor.write('SOUR:AM:INT:FREQ {}'.format(freq))
        # self._thor.write('SOUR:AM:INT[:DEPT] <{}>'.format(depth)) # NO FUNCIONA
        
    def offModulation(self): # no anda
        self._thor.write('SOUR:AM[:STAT] OFF')
