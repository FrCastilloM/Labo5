# -*- coding: utf-8 -*-
"""
Created on Feb 2024

@author: Francisco Castillo - Mauro Chavez - Tomas Spinazzola

"""
import matplotlib.pyplot as plt
import pyvisa as visa
import numpy as np
import time
import re
from operator import itemgetter

import os
os.chdir("Z:\\Fisica\\Exactas\\Laboratorio 5\\Datos y Analisis\\1. Glow")

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt5')

plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)

from TektronixAFG3021B import AFG3021B

#%% Comunicacion con instrumentos

rm = visa.ResourceManager()

rm.list_resources()

mult1 = rm.open_resource('GPIB0::22::INSTR') 
mult2 = rm.open_resource('GPIB0::23::INSTR') 
fungen = AFG3021B('USB0::0x0699::0x0346::C036492::INSTR')

#%% #%% Barrido en voltajes

def graficar_barrido(j = 1, xscale = 'linear'):
    
    voltajes, v1, v2 = np.loadtxt('barrido_voltajes_{}.txt'.format(j), delimiter=',', unpack=True, skiprows = 1)

    v_ac = v1*1023.1 - v2   # [V]
    I = v2/14.86 # [mA]
    n = len(I)//2

    plt.figure(j)
    plt.title('Curva V vs I',fontsize=14)
    plt.xscale(xscale)
    plt.plot(I[:n], v_ac[:n], '.', color="tab:blue", label = 'Subida')
    plt.plot(I[n:], v_ac[n:], 'x', color="tab:red", label = 'Bajada')
    plt.xlabel('Corriente [mA]',fontsize=14)
    plt.ylabel('Voltaje anodo-catodo [V]',fontsize=14)
    plt.legend(fontsize=14, loc = 'lower right')
    plt.tight_layout()
    plt.grid(True)

def barrido_voltajes(j = 1, n = 10, v_i = 0, v_f = 800, fc = 500, ts = 1, xscale = 'linear', p = 'None', d = 'None'):
    
    # j = etiqueta de la medicion
    # n = numero de puntos (subida/bajada)
    # v_i = voltaje inicial (fuente alta tension)
    # v_f = voltaje final (fuente alta tension)
    # fc = factor de calibracion generador-fuente
    # ts = tiempo de espera entre seteo de voltaje y medicion
    # p = presion
    # d = gap interelectrodico
    
    name = 'barrido_voltajes_{}.txt'.format(j)
    
    if os.path.isfile(name):
        print('Error: cambiar j.')
        return

    # barrido de voltajes subida y bajada
    voltajes = np.hstack((np.linspace(v_i, v_f, n)/fc, np.linspace(v_f, v_i, n)/fc))     

    v1 = np.zeros(len(voltajes)) # mult1 = "tension" 
    v2 = np.zeros(len(voltajes)) # mult2 = "corriente"

    for i in range(len(voltajes)):
        fungen.setOffset(voltajes[i], channel = 1)
        time.sleep(ts)
        v1[i], v2[i] = float(mult1.query('MEASURE:VOLTAGE:DC?')), float(mult2.query('MEASURE:VOLTAGE:DC?'))
        print(i)
        
    fungen.setOffset(0, channel = 1)

    
    np.savetxt(name, np.column_stack([voltajes, v1, v2])
                ,delimiter=',', header = 'Vin [V], V1 [V], V2 [V], p = {pres} Torr, d = {dist} cm, ts = {ts} s, n = {n}, fc = {fc}, vi = {vi}, vf = {vf}'.format(pres = p, dist = d, ts = ts, n = n, fc = fc, vi = v_i, vf = v_f))
    
    graficar_barrido(j, xscale)

#%% #%% Voltaje de ruptura

def graficar_vruptura(j = 1, xscale = 'linear'):
    
    voltajes, v1, v2 = np.loadtxt('paschen_{}.txt'.format(j), delimiter=',', unpack=True, skiprows = 1)

    v_ac = v1*1023.1 - v2   # [V]
    I = v2/14.86 # [mA]

    V_I = np.zeros((len(v_ac),2))
    V_I[:,0] = v_ac
    V_I[:,1] = I
    vc = max(filter(lambda a: a[1] <= 1e-3, V_I), key=itemgetter(0))

    plt.figure(j)
    plt.title('Curva V vs I',fontsize=14)
    plt.xscale(xscale)
    plt.scatter(vc[1], vc[0], c = 'tab:red')
    plt.plot(I, v_ac, '.', color="tab:blue")
    plt.xlabel('Corriente [mA]',fontsize=14)
    plt.ylabel('Voltaje anodo-catodo [V]',fontsize=14)
    plt.legend(fontsize=14)
    plt.tight_layout()
    plt.grid(True)


'''

El siguiente codigo ejecuta, idealmente, dos barridos, el primero es de testeo para ubicar la region
del voltaje de ruptura, el segundo es un barrido mas fino para medirlo con una buena presicion

'''

def vruptura(j = 1, n_a = 15, n_b = 40, vi = 200, vf = 900, fc = 500, ts = 1, xscale = 'linear', p = 'None', d = 'None'):
    
    # j = etiqueta de la medicion
    # n_a = numero de puntos barrido de testeo (solo subida)
    # n_b = numero de puntos barrido denso
    # vi = voltaje inicial (fuente alta tension)
    # vf = voltaje final (fuente alta tension)
    # fc = factor de calibracion generador-fuente
    # ts = tiempo de espera entre seteo de voltaje y medicion
    # p = presion
    # d = gap interelectrodico
    
    name = 'paschen_{}.txt'.format(j)
    
    if os.path.isfile(name):
        print('Error: cambiar j.')
        return

    # barrido de testeo (subida)
    voltajes = np.linspace(vi,vf,n_a)/fc    

    R = 0
    while R < 2:
        v1 = np.zeros(len(voltajes)) # divisor tension 
        v2 = np.zeros(len(voltajes)) # caida en R de corriente
        
        fungen.setOffset(0, channel = 1)
        time.sleep(1)
        
        for i in range(len(voltajes)):
            fungen.setOffset(voltajes[i], channel = 1)
            time.sleep(ts)
            v1[i], v2[i] = float(mult1.query('MEASURE:VOLTAGE:DC?')), float(mult2.query('MEASURE:VOLTAGE:DC?'))
            print('j = ', i)
            
            v_ac = v1[i]*1023.1 - v2[i]
            I_i = v2[i]/14.86
            print('vac = ', v_ac, 'V')
            print('i = ', I_i, 'mA')
            print('vfuente = ', v1[i]*1023.1, 'V')
            print(' ')
            if I_i > 1e-3:
                if R == 0:
                    voltajes = np.hstack((np.linspace(voltajes[i-4],voltajes[i-2],5),
                                          np.linspace(voltajes[i-2],voltajes[i+1],n_b)))
                print('--break--')
                print(' ')
                break
        R += 1
        
    fungen.setOffset(0, channel = 1)

    np.savetxt('paschen_aire_{}.txt'.format(j),
                np.column_stack([voltajes, v1, v2])
                ,delimiter=',', header = 'Vin [V], V1 [V], V2 [V], p = {pres} Torr, d = {dist} cm, ts = {tss} s, n_a = {na}, n_b = {nb}, fc = {fcc}, vi = {vii}, vf = {vff}'.format(pres = p, dist = d, tss = ts, na = n_a, nb = n_b, fcc = fc, vii = vi, vff = vf))

    graficar_vruptura(j, xscale)
    


#%% Curva de Pashen

# mediciones descartadas
med_throw = []

v_ruptura = []
pd = []

for l in range(1,22):
    if l in med_throw:
        continue
    archivo = 'paschen_{}.txt'.format(l)
    vi, v1, v2 = np.loadtxt(archivo, delimiter=',', unpack=True, skiprows = 1)
    with open(archivo, 'r') as file:
        header = file.readline() 
    h = re.findall(r'\d+\.\d+|\d+', header)
    p = float(h[2]) # Torr
    d = float(h[3]) # cm
        
    v_ac = v1*1023.1 - v2
    I = v2/14.86

    V_I = np.zeros((len(v_ac),2))
    V_I[:,0] = v_ac
    V_I[:,1] = I
    print(f'j={l}')
    if all(I <= 1e-10) or all(I >= 1e-10): # Para que saltee los casos en los que no llego a saturar
        continue
    vc = max(filter(lambda a: a[1] <= 1e-3, V_I), key=itemgetter(0))
    print(vc)
    v_ruptura.append([vc[0]])
    
    # plt.figure(l)
    # plt.title('Curva V vs I',fontsize=14)
    # # plt.xscale('linear')
    # plt.xscale('log')
    # plt.scatter(vc[1], vc[0], c = 'tab:red')
    # plt.plot(I, v_ac, '.', color="tab:blue")
    # plt.xlabel('Corriente [mA]',fontsize=14)
    # plt.ylabel('Voltaje anodo-catodo [V]',fontsize=14)
    # plt.legend(fontsize=14)
    # plt.tight_layout()
    # plt.grid(True)
    
    pd.append(d*p)
    
plt.figure(99)
plt.title('Curva de Pashen',fontsize=14)
plt.plot(pd, v_ruptura, '.', label = 'Aire')
plt.xlabel(r'pd [cm$\cdot$Torr]',fontsize=14)
plt.ylabel('Voltaje ruptura [V]',fontsize=14)
plt.legend(fontsize=14)
plt.tight_layout()
plt.grid(True)

