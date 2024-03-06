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
daq_ch2 = 'Dev2/ai2'

num_ai_channels = 2
ai_channels = [1, 2]

# Configuración de la tarea
    
task = nidaqmx.Task()
# task.close()

# Setear el modo y rango de un canal analógico
modo = nidaqmx.constants.TerminalConfiguration(10083)
ai_ch1 = task.ai_channels.add_ai_voltage_chan(daq_ch1, terminal_config = modo, max_val = 10, min_val = -10)
ai_ch2 = task.ai_channels.add_ai_voltage_chan(daq_ch2, terminal_config = modo, max_val = 10, min_val = -10)


print('')
print("Número de canales AI configurados:", num_ai_channels)
print("Canales AI configurados:", ai_channels)
print('')
print('------ Channel 1 ------')
print(ai_ch1.ai_term_cfg)
print(ai_ch1.ai_max)
print(ai_ch1.ai_min)
print('------ Channel 2 ------')
print(ai_ch2.ai_term_cfg)
print(ai_ch2.ai_max)
print(ai_ch2.ai_min)
print('-----------------------')

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

''' espectro_gamma_2(name, dt, cant_ventanas, task)

Extension de la funcion espectro_gamma de adquisicion_nuclear.py para medir con 
un numero arbitrario de canales (>1). 

Los parametros umbral, HV, distance y sign ahora son listas donde cada posicion
corresponde a un canal del daq diferente, en el orden en el que fueron añadidos

num_canales = numero de canales conectados

Advertencia!! dado que la frecuencia maxima del daq es 4e5 kS/s, este valor 
se tiene que dividir entre el numero de canales utilizados, es decir, si hay
dos canales conectados, para mantener la cantidad maxima de data x seg. que se puede tomar
fs tiene que ser fs = 2e5

'''


def espectro_gamma_2(name, dt, cant_ventanas, umbral, task, HV = ['None', 'None'],
                  save = 'False', fs = 2e5 , distance = [10, 10],
                   sign = [-1, 1], remplace = 'False', num_canales = 2):

    # Verifico que no haya guardado un archivo con el mismo nombre
    if os.path.isfile(name) and remplace != 'True':
        print('Error: cambiar etiqueta o indicar remplace = True.')
        return
    
    # Chequeo que las dimensiones de todos los parametros coincidan
    if any(x != num_canales for x in [len(umbral), len(distance), len(sign), len(HV)]):
        print('Error: mismatch entre dimensiones de los parametros con el numero de canales.')
        
    # Lista de listas [[ch1], [ch2], ...]
    picos = [list() for _ in range(num_canales)] # Tensiones de los picos
    tiempos = [list() for _ in range(num_canales)] # Posiciones temporales de los picos
    dist = [list() for _ in range(num_canales)] # Distancia temporal entre los picos
    eventos_por_ventana = [list() for _ in range(num_canales)] # Numero de eventos por ventana
    
    print('----- Espectro gamma: {} -----'.format(name))
    for i in range(cant_ventanas):
    
        # Medimos en la ventana temporal
        y = np.diag(sign) @ medicion_una_vez(dt, fs, task)
        t = np.arange(len(y[0,:]))/fs
        
        # (opcional) Guardo datos de la ventana temporal
        if save == 'True':
            np.savetxt(name + '_{}.txt'.format(i), np.column_stack([list(t)] + y.tolist())
                        ,delimiter=',',
                        header = 'Tiempo [s], Voltaje [V] ..., frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV))
            
        print(i+1, 'de', cant_ventanas, end = '')
        for q in range(num_canales):
            
            # Deteccion de picos, umbral y distancia son parametros opcionales de la funcion
            peaks, prop = find_peaks(y[q,:], 
                                     height = umbral[q], width = 0, distance = distance[q])        
        
            # Agrego datos a las listas
            picos[q].extend(y[q,:][peaks])
            tiempos[q].extend(t[peaks])
            dist[q].extend(np.diff(peaks)/fs)
            eventos_por_ventana[q].extend([len(peaks)])

            print('       # eventos (ai{})= '.format(ai_channels[q]), len(peaks), end = '')
        print('')
    
    # Guardo resultados finales
    for q in range(num_canales):
        np.savetxt(name + '_ai{}_ '.format(ai_channels[q]) + '_resultados.txt', 
                    np.column_stack([tiempos[q], picos[q]]), delimiter=',',
                    header = 'Tiempo(peaks) [s], Voltaje(peaks) [V], frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV[q]))
        np.savetxt(name + '_ai{}_ '.format(ai_channels[q]) + '_resultados_dist.txt', 
                   dist[q], header = 'Distancia temporal entre picos [s]')
        np.savetxt(name + '_ai{}_ '.format(ai_channels[q]) + '_resultados_eventos.txt', 
                   eventos_por_ventana[q], header = 'Numero de eventos detectados en ventanas de dwell time = {dt} s'.format(dt = dt))
        
    print('')
    print('  /\_/\ ')
    print(' ( o.o )')
    print('  > ^ < ')
    print('----------------------------------------------')
    
    return picos, tiempos, dist, eventos_por_ventana


''' coincidencias(name, dt, cant_ventanas, tau_criterio, alpha, task)

Esta funcion toma la data del daq por ventanas para dos canales, busca los picos y analiza si
hay algun pico del canal 1 que este cerca temporalmente de un pico del canal 2. Esto se toma
como una coincidencia. El codigo cuenta y guarda el numero de coincidencias por ventana, la 
amplitud de los picos de coincidencia (de cada canal), y los tiempos.

Por otro lado, tambien se toman y guardan los datos el espectro gamma para c/canal.

tau_criterio = picos separados temporalmente por debajo de este valor se cuentan como una coincidencia

alpha = angulo relativo a la fuente del centellador desplazable (este valor se anota en el header
de los .txt con los datos)
                                                          
fotopico_volt = [[min_ai1, max_ai1], [min_ai2, max_ai2]]
Si la amplitud en voltaje de un pico detectado en el canal ai1(2) esta dentro del rango establecido
por este parametro, entonces es candidato a ser considerado como una coincidencia. De lo contrario,
si el pico cae por fuera de este rango, es descartado.
Estos rangos deben estar centrados en los fotopicos medidos para cada canal correspondientes a la emision
de 0.511 MeV del 22Na.   

'''

def criterio_temporal(picoai1, picoai2, fs, tau_criterio):
    return abs((picoai1 - picoai2)/fs) < tau_criterio

def criterio_voltaje(picoai1_v, picoai2_v, umbrales):
    return umbrales[0][0] <= picoai1_v <= umbrales[0][1] and umbrales[1][0]  <= picoai2_v <= umbrales[1][1] 


def coincidencias(name, dt, cant_ventanas, tau_criterio, alpha, fotopico_volt, task, HV = ['None', 'None'],
                  save = 'False', fs = 2e5, umbral = [[0.05,9.5], [0.05,9.5]] , distance = [10, 10],
                   sign = [-1, 1], remplace = 'False', num_canales = 2):

    # Verifico que no haya guardado un archivo con el mismo nombre
    if os.path.isfile(name) and remplace != 'True':
        print('Error: cambiar etiqueta o indicar remplace = True.')
        return
    
    # Chequeo que las dimensiones de todos los parametros coincidan
    if any(x != num_canales for x in [len(umbral), len(distance), len(sign), len(HV)]):
        print('Error: mismatch entre dimensiones de los parametros con el numero de canales.')
        
    # Lista de listas [[ch1], [ch2], ...]]
    picos = [list() for _ in range(num_canales)] # Tensiones de los picos
    tiempos = [list() for _ in range(num_canales)] # Posiciones temporales de los picos
    dist = [list() for _ in range(num_canales)] # Distancia temporal entre los picos
    eventos_por_ventana = [list() for _ in range(num_canales)] # Numero de eventos por ventana
    
    coin_picos = [list() for _ in range(num_canales)] # Tensiones de los picos
    coin_tiempos = [list() for _ in range(num_canales)] # Posiciones temporales de los picos
    coin_dist = [list() for _ in range(num_canales)] # Distancia temporal entre los picos
    coincidencias_por_ventana = []
    
    print('----- Espectro gamma - coincidencias: {} -----'.format(name))
    for i in range(cant_ventanas):
    
        # Medimos en la ventana temporal
        y = np.diag(sign) @ medicion_una_vez(dt, fs, task)
        t = np.arange(len(y[0,:]))/fs
        
        # (opcional) Guardo datos de la ventana temporal
        if save == 'True':
            np.savetxt(name + '_{}.txt'.format(i), np.column_stack([list(t)] + y.tolist())
                        ,delimiter=',',
                        header = 'Tiempo [s], Voltaje [V] ..., frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV))
        
        picos_temp = [list() for _ in range(num_canales)]
        print(i+1, 'de', cant_ventanas, end = '')
        for q in range(num_canales):
            
            # Deteccion de picos, umbral y distancia son parametros opcionales de la funcion
            peaks, prop = find_peaks(y[q,:], 
                                     height = umbral[q], width = 0, distance = distance[q])        
            
            picos_temp[q].extend(peaks)
            
            # Agrego datos a las listas
            picos[q].extend(y[q,:][peaks])
            tiempos[q].extend(t[peaks])
            dist[q].extend(np.diff(peaks)/fs)
            eventos_por_ventana[q].extend([len(peaks)])

            print('  # event. (ai{}) = '.format(ai_channels[q]), len(peaks), end = '')
        
        # Cuento las coincidencias
        coincidences_index = [list() for _ in range(num_canales)]
        
        for v in range(len(picos_temp[0])):
            for mu in range(len(picos_temp[1])):
                ct = criterio_temporal(picos_temp[0][v], picos_temp[1][mu], fs, tau_criterio)
                cv = criterio_voltaje(y[0,:][picos_temp[0][v]], y[1,:][picos_temp[1][mu]], fotopico_volt)
   
                if ct and cv:
                    coincidences_index[0].extend([picos_temp[0][v]])
                    coincidences_index[1].extend([picos_temp[1][mu]])
                    
        for q in range(num_canales):
            coin_picos[q].extend(y[q,:][coincidences_index[q]])
            coin_tiempos[q].extend(t[coincidences_index[q]])
            coin_dist[q].extend(np.diff(coincidences_index[q])/fs)
        
        coincidencias_por_ventana.append(len(coincidences_index[0]))
        
        print('  # coinci. = ', len(coincidences_index[0]))
        
    
    # Guardo resultados finales
    for q in range(num_canales):
        np.savetxt(name + '_ai{}'.format(ai_channels[q]) + '_resultados.txt', 
                    np.column_stack([tiempos[q], picos[q]]), delimiter=',',
                    header = 'Tiempo(peaks) [s], Voltaje(peaks) [V], frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV[q]))
        np.savetxt(name + '_ai{}'.format(ai_channels[q]) + '_resultados_dist.txt', 
                   dist[q], header = 'Distancia temporal entre picos [s]')
        np.savetxt(name + '_ai{}'.format(ai_channels[q]) + '_resultados_eventos.txt', 
                   eventos_por_ventana[q], header = 'Numero de eventos detectados en ventanas de dwell time = {dt} s'.format(dt = dt))
    
    np.savetxt(name + '_coincidencias_resultados.txt', 
                np.column_stack([coin_tiempos[0], coin_picos[0], coin_tiempos[1], coin_picos[1]]), delimiter=',',
                header = 'Tiempo(coincidence-peaks-ch1) [s], Voltaje(coincidence-peaks-ch1) [V], Tiempo(coincidence-peaks-ch2) [s], Voltaje(coincidence-peaks-ch2) [V], frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV, angulo = {alph}'.format(fs = fs, dt = dt, HV = HV, alph = alpha))
    
    np.savetxt(name + '_coincidencias_resultados_dist.txt', 
               np.column_stack([coin_dist[0], coin_dist[1]]), header = 'Distancia temporal entre picos: ch1 [s], ch2 [s], angulo = {alph}'. format(alph = alpha))
    
    np.savetxt(name + '_coincidencias_resultados_eventos.txt', 
               coincidencias_por_ventana, header = 'Numero de coincidencias en ventanas de dwell time = {dt} s, angulo = {alph}'.format(dt = dt, alph = alpha))

    
    print('')
    print('  /\_/\ ')
    print(' ( o.o )')
    print('  > ^ < ')
    print('----------------------------------------------')
    
    return picos, tiempos, dist, eventos_por_ventana, coin_picos, coin_tiempos, coin_dist, coincidencias_por_ventana

''' graficar(picos, dist, eventos, channel, dt)

Funcion simple para graficar. Esta funcion es mas que nada para reducir un poco
de espacio en el codigo.

Picos, dist y eventos son la data que sale de espectro_gamma
dt = duracion de la medicion
channel = canal de la data

'''

def graficar(picos, dist, eventos, channel, dt):
    
    plt.figure()
    plt.title('Amplitud de los picos', fontsize = 14)
    plt.hist(picos, bins=300, label = 'ai{}'.format(channel))
    plt.xlabel('Tensión [V]', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)
    
    plt.figure()
    plt.title('Distancia entre picos', fontsize = 14)
    plt.hist(dist, bins=200, label = 'ai{}'.format(channel))
    plt.xlabel(r'$\Delta$t [s]', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)
    
    plt.figure()
    plt.title('Detecciones en dt = {} s'.format(dt), fontsize = 14)
    plt.hist(eventos, bins='auto', label = 'ai{}'.format(channel))
    plt.xlabel('n', fontsize = 14)
    plt.ylabel('Cuentas', fontsize = 14)
    plt.yscale("log")
    plt.legend(fontsize = 14)

#%% Medicion simple de testeo

HV = [1.02, 1.02]
fs = 2e5
dt = 5
umbral = [[.05, 9.5], [.05, 9.5]]
tau_criterio = 60*1e-6 # s

y = np.diag([-1,1]) @ medicion_una_vez(dt, fs, task)
t = np.arange(len(y[0,:]))/fs

# np.savetxt('data_daq_prueba_2.txt', np.column_stack([list(t)] + y.tolist())
#             ,delimiter=',',
#             header = 'Tiempo [s], Voltaje [V] ..., frecuencia de muestreo = {fs} Hz, dwell time = {dt} s, HV = {HV} kV'.format(fs = fs, dt = dt, HV = HV))

plt.figure()
plt.title("Señal DAQ",fontsize=14)
plt.ylabel(r'Voltaje [V]',fontsize=14)
plt.xlabel(r'Tiempo [s]',fontsize=14)

picos_temp = [list() for _ in range(2)]
for q in range(num_ai_channels):
    plt.plot(t, y[q,:], label = 'ai{}'.format(ai_channels[q]))
    picos, _ = find_peaks(y[q,:], height = umbral[q], width = 0, distance = 10)   
    picos_temp[q].extend(picos)     
    plt.plot(t[picos],y[q,:][picos], 'x', label = 'peaks')
    print(len(picos))

coincidences_index = [list() for _ in range(2)]

for v in range(len(picos_temp[0])):
    for mu in range(len(picos_temp[1])):
        if abs((picos_temp[0][v] - picos_temp[1][mu])/fs) < tau_criterio:
            coincidences_index[0].extend([picos_temp[0][v]])
            coincidences_index[1].extend([picos_temp[1][mu]])
            
print('coincidencias: ', len(coincidences_index[0]))

for i in range(len(coincidences_index[0])):
    plt.axvline(x = t[coincidences_index[0][i]], linestyle = '--', c = 'red', lw = 1)
    plt.axvline(x = t[coincidences_index[1][i]], linestyle = '--', c = 'k', lw = 1)
    
plt.legend(fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()


#%% Medicion definitiva, construccion del histograma

HV = [1.02, 1.02] # kV
dt = .5 # dwell time/duracion
n = 1000 # cantidad de ventanas 
umbral = [[.05, 9.5], [.05, 9.5]] # umbrales de deteccion de picos

name = '207Cs'

picos, tiempos, dist, eventos_por_ventana = espectro_gamma_2(name, dt, n, umbral, task, HV = HV, 
                                                           save = 'False', sign = [-1, 1]) 

for q in range(num_ai_channels):
    graficar(picos[q], dist[q], eventos_por_ventana[q], ai_channels[q], dt)


#%% Cargar datos y graficar

name = '207Cs_1_ai1_resultados'

tiempos, picos = np.loadtxt(name + '.txt', delimiter=',', unpack=True, skiprows = 1)
dist = np.loadtxt(name + '_dist.txt', skiprows = 1)
eventos = np.loadtxt(name + '_eventos.txt', skiprows = 1)

graficar(picos, dist, eventos, 1, .5) 

#%% Coincidencias

alpha = 33 # grados
tau_criterio = 60*1e-6 # s
HV = [1.02, 1.02] # kV
dt = 5 # dwell time/duracion
n = 100 # cantidad de ventanas 
umbral = [[.05, 9.5], [.05, 9.5]] # umbrales de deteccion de picos

# regiones de las campanas definidas por los fotopicos de 0.511 MeV 
fotopico_volt = [[0.86, 1.5], [1.8, 2.8]] 

name = '22Na_coincidencias'

(picos, tiempos, dist, eventos_por_ventana, coin_picos, coin_tiempos, 
    coin_dist, coincidencias_por_ventana) = coincidencias(name, dt, n, tau_criterio, 
                                                          alpha, fotopico_volt, task, HV = HV,
                                                          save = 'False', umbral = umbral)

                                                          
for q in range(num_ai_channels):
    graficar(picos[q], dist[q], eventos_por_ventana[q], ai_channels[q], dt)    

q = 0
graficar(coin_picos[q], coin_dist[q], coincidencias_por_ventana, 
         str(ai_channels[q]) + ' coincidences', dt)





