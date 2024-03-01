# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 21:55:55 2024

@author: Francisco
"""

import matplotlib.pyplot as plt
import numpy as np
import nidaqmx
from scipy.signal import find_peaks

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt5')

import os
# reemplazar por la carpeta donde quiere guardar los datos
os.chdir("C:\\Francisco\\Cursando notebook\\L5\\nuclear\\dia4") 

plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)

#%% Mediciones con DAQ

# para saber el ID de la placa conectada (DevX)
system = nidaqmx.system.System.local()
for device in system.devices:
    print(device)

daq_ch1 = 'Dev2/ai1' 

# Configuración de la tarea
    
task = nidaqmx.Task()

# Setear el modo y rango de un canal analógico
modo = nidaqmx.constants.TerminalConfiguration(10083)
ai_ch1 = task.ai_channels.add_ai_voltage_chan(daq_ch1, terminal_config = modo, max_val = 10, min_val = -10)


print('------ Channel 1 ------')
print(ai_ch1.ai_term_cfg)
print(ai_ch1.ai_max)
print(ai_ch1.ai_min)
print('-----------------------')

# por si alguna vez el daq falla y necesitan cerrar la tarea y volver a iniciarla:
# task.close()

#%% Funciones para medir

''' medicion_una_vez(duracion, fs, task)

Medicion simple para obtener los datos adquiridos por la daq

parametros:
    
duracion = duracion de la medicion

fs = frecuencia de muestreo (usar siempre la max permitida 4e5)

task = tarea del daq iniciada (ver celda anterior)

'''

def medicion_una_vez(duracion, fs, task):
    
    # Cantidad de puntos de la medición
    cant_puntos = int(duracion*fs)
        
    # Frecuencia de sampleo:
    task.timing.cfg_samp_clk_timing(rate = fs, samps_per_chan = cant_puntos,
                                    sample_mode = nidaqmx.constants.AcquisitionType.FINITE)
        
    # Adquirimos datos:
    datos = task.read(number_of_samples_per_channel = nidaqmx.constants.READ_ALL_AVAILABLE, 
                      timeout = nidaqmx.constants.WAIT_INFINITELY) 
    # Sin el timeout no podemos cambiar la duracion a mas de 10 seg.
        
    datos = np.asarray(datos)    
    return datos # Tenemos un matriz de [# de canales]x[cant_puntos], donde
                    # cada fila son los datos de todo un canal.


''' espectro_gamma(name, dt, cant_ventanas, task)

Esta funcion es para adquirir el espectro gamma emitido por una fuente,
el codigo toma tiras de datos en la daq y anota la altura de los pulsos, su
posicion temporal, la distancia relativa entre ellos, y el numero de pulsos por 
ventana temporal 


parametros importantes:

name = nombre del archivo .txt que se va a guardar con los datos

dt = duracion de la medicion (parametro duracion de la funcion anterior)

cantidad_ventanas = numero de veces que el daq va a tomar una tira de datos
recomendamos usar cant_ventanas = 3000 junto con dt = 0.5 para tener una buena estadistica (25m)

umbral = parametro height de find_peaks = umbral donde el codigo buscará picos
recomendamos tomar datos crudos con el daq y definir un valor de umbral por encima
del ruido electrico y por debajo del limite del daq (10 V)

parametros opcionales:
    
HV = High Voltaje, por si quieren guardar en el header del .txt con los datos 
el valor de voltaje que tienen en la fuente

save = condicion por si quieren guardar todas las ventanas de datos (voltaje vs tiempo),
advertencia! si la frecuencia es de 4e5 y la duracion de .5, cada .txt ocupa unos 10MB de 
memoria, ademas de volver mas lento el codigo. Recomendamos usar esta opcion solo para
guardar algunos pocos datos.

distance = parametro distance de la funcion height de find_peaks = distancia entre datos
minima para deteccion de un pico

sign = +/- 1 que multiplica a los datos
si se esta usando un centellador de polaridad negativa poner -1, de lo contrario 1

remplace = por si ya existe un .txt con el name indicado y se lo quiere reemplazar

'''

def espectro_gamma(name, dt, cant_ventanas, umbral, 
                   task, fs = 4e5, HV = 'None', save = 'False', distance = 10,
                   sign = -1, remplace = 'False'):
    
    # Verifico que no haya guardado un archivo con el mismo nombre
    if os.path.isfile(name) and remplace != 'True':
        print('Error: cambiar etiqueta o indicar remplace = True.')
        return
    
    picos = np.zeros(0) # Tensión del pico
    tiempos = np.zeros(0) # Posición temporal del pico
    dist = np.zeros(0) # Distancia temporal entre los picos
    eventos_por_ventana = [] # Numero de eventos por venana
    
    print('----- Espectro gamma: {} -----'.format(name))
    for i in range(cant_ventanas):
    
        # Medimos en la ventana temporal
        y = sign * medicion_una_vez(dt, fs, task)
        t = np.arange(len(y))/fs
        
        # (opcional) Guardo datos de la ventana temporal
        if save == 'True':
            np.savetxt(name + '_{}.txt'.format(i), np.column_stack([t, y])
                        ,delimiter=',',
                        header = 'Tiempo [s], Voltaje [V], frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV))
            
        # Deteccion de picos. Umbral y distancia son parametros opcionales de la funcion
        peaks, prop = find_peaks(y, height = umbral, width = 0, distance = distance)        
        
        # Agrego datos a las listas
        picos = np.concatenate((picos, y[peaks]))
        tiempos = np.concatenate((tiempos, t[peaks]))
        dist = np.concatenate((dist, np.diff(peaks)/fs))
        eventos_por_ventana.append(len(peaks))

        print(i+1, "de", cant_ventanas,"       # eventos = ", len(peaks))
    
    # Guardo resultados finales
    np.savetxt(name + '_resultados.txt', np.column_stack([tiempos, picos])
                ,delimiter=',',
                header = 'Tiempo(peaks) [s], Voltaje(peaks) [V], frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV))

    np.savetxt(name + '_resultados_dist.txt', dist, header = 'Distancia temporal entre picos [s]')
    np.savetxt(name + '_resultados_eventos.txt', eventos_por_ventana, header = 'Numero de eventos detectados en ventanas de dwell time = {dt} s'.format(dt = dt))
    
    print('')
    print("  /\_/\ ")
    print(" ( o.o )")
    print("  > ^ < ")
    print('----------------------------------------------')
    
    return picos, tiempos, dist, eventos_por_ventana

''' graficar(picos, dist, eventos, dt)

Funcion simple para graficar. Esta funcion es mas que nada para reducir un poco
de espacio en el codigo.

Picos, dist y eventos son la data que sale de espectro_gamma
dt = duracion de la medicion

'''

def graficar(picos, dist, eventos, dt):
    
    plt.figure()
    plt.title('Amplitud de los picos', fontsize = 14)
    plt.hist(picos, bins=300, label = 'data')
    plt.xlabel('Tensión [V]', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)
    
    plt.figure()
    plt.title('Distancia entre picos', fontsize = 14)
    plt.hist(dist, bins=200, label = 'data')
    plt.xlabel(r'$\Delta$t [s]', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)
    
    plt.figure()
    plt.title('Detecciones en dt = {} s'.format(dt), fontsize = 14)
    plt.hist(eventos, bins='auto', label = 'data')
    plt.xlabel('n', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)
    
#%% Medicion simple de testeo

fs = 4e5 # Frecuencia de muestreo
dt = .5 # segundos de medicion

y = -medicion_una_vez(dt,fs,task)
t = np.arange(len(y))/fs

# t, y = np.loadtxt('207Bi_1_0.txt', delimiter=',', unpack=True, skiprows = 1)

picos, _ = find_peaks(y, height = [0.05, 9.5], width = 0, distance = 10)

plt.figure()
plt.plot(t,y,label="Datos")
plt.title("Señal DAQ",fontsize=14)
plt.ylabel(r'Voltaje [V]',fontsize=14)
plt.xlabel(r'Tiempo [s]',fontsize=14)
plt.plot(t[picos],y[picos], 'x', label = 'peaks')
plt.legend(fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()

print(len(picos))

#%% Medicion definitiva, construccion del histograma

HV = 1.02 # kV
dt = .5 # dwell time/duracion
n = 3000 # cantidad de ventanas 
umbral = [0.05, 9.5] # [cota inf, cota sup]

name = '133Ba'

picos, tiempos, dist, eventos_por_ventana = espectro_gamma(name, dt, n, umbral, task, HV = HV,
                                                           save = 'False', sign = -1) 

graficar(picos, dist, eventos_por_ventana, dt)

#%% Cargo datos y grafico

name = '22Na_1'

tiempos, picos = np.loadtxt(name + '_resultados.txt', delimiter = ',', unpack = True, skiprows = 1)
dist = np.loadtxt(name + '_resultados_dist.txt', skiprows = 1)
eventos = np.loadtxt(name + '_resultados_eventos.txt', skiprows = 1)





