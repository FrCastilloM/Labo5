# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 11:09:40 2024

@author: Francisco
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'qt5')

import os

PC = 1 

if PC==1:
    os.chdir("Z:\\Fisica\\Codigos")
    from DataAnalisys import stats
    os.chdir("Z:\\Fisica\\Exactas\\Laboratorio 5\\Datos y Analisis\\2. Nuclear\\dia3")
if PC==2:
    os.chdir("C:\\Francisco\\Cursando notebook\\L5\\nuclear\\dia3")

plt.rc('xtick', labelsize=13)
plt.rc('ytick', labelsize=13)

#%% Funciones

def contar_eventos(times, dt, dt_old):    
    eventos = []
    a = int(dt_old/dt)
    i = 0
    while 1:
        R = np.zeros(a)
        for j in range(a):
            while dt*j <= times[i] < dt*(j+1):
                i += 1
                if i == len(times):
                    R[j] = i
                    eventos.extend(R)
                    return np.hstack((eventos[0], np.diff(eventos)))
            R[j] = i
        eventos.extend(R)

def stirling(n):
    if n>10:
        return n*np.log(n)-n+.5*np.log(2*np.pi*n)+ 1/(12*n)
    else:
        return np.log(np.math.factorial(n))

def Poisson(k, mu):
    P = np.zeros(len(k))
    for i in range(len(k)):
       P[i] = mu**int(k[i])*np.exp(-mu)/np.exp(stirling(int(k[i])))
    return P

# defino la distribucion de Polya-Aeppli 
def PA(x,mu,p):
    P = np.zeros(len(x))
    for i in range(len(x)):
        n = int(x[i])
        if n == 0:
            P[i] = np.exp(-mu)
        if n>0:
            suma = 0
            for k in range(1, n+1):
                suma += mu**k * (1-p)**(n-k) * p**k * np.exp(stirling(n-1)-stirling(n-k)-stirling(k-1)-stirling(k))
            P[i] = np.exp(-mu)*suma
    return P

def contar_repeticiones(lista, numero):
    contador = 0
    for i in range(len(lista)):
        if lista[i] == numero:
            contador += 1
    return contador

# estimador de la desviacion estandar
def sdev(X):
    return np.sqrt(np.sum((X-np.mean(X))**2)/(len(X)-1))

def Hist(n, N, K, g):
    R = np.zeros((K,n+1)) # genero una matriz K x n+1 que contendra todos los resultados
    # completo la matriz R
    for k in range(K):
        reps = np.zeros(n+1) 
        for j in range(n+1):
            reps[j] = contar_repeticiones(g[N*k:N*(k+1)],j) # cuento las veces que se obtuvo un numero j de exitos en los N experimentos
        R[k,:] = reps 
        
    # la matriz R asi construida tiene en sus filas -> R[k,:] = (# veces que detecte 0 fotones, # veces que detecte 1 foton, ...)
    # habiendo repetido el experimento N veces. Luego, cada una de sus columnas es una repeticion de esto ultimo.
    
    # ahora usamos las diferentes filas de R para obtener el valor medio X, y su desviacion estandar errX
    X = np.zeros(n+1)
    errX = np.zeros(n+1)
    # calculo los estimadores para cada bin
    for i in range(n+1):
        X[i] = np.mean(R[:,i])
        errX[i] = 2*sdev(R[:,i])/np.sqrt(len(R[:,i])) # error del promedio (2 sigmas)
    return X, errX

def Ajuste_poisson(bin_centers, X, N, errX, p0):
    bc = bin_centers[X>0]
    parameters, covariance = curve_fit(Poisson, bc, X[bc]/N, sigma = errX[bc]/N, 
                                       p0 = [p0], absolute_sigma = True)
    paramErrors = np.sqrt(np.diag(covariance))
    
    print("-------Ajuste Poisson------")
    
    parametros=['mu']
    for xe in range(len(parametros)):
        print(parametros[xe],"=",parameters[xe],"+/-",paramErrors[xe])
        
    mu = parameters 
    
    ajuste = Poisson(bc, mu)
    
    EstadisticaDelAjuste = stats(ajuste, bc, X[bc]/N, 1)
    
    print("R^2:",EstadisticaDelAjuste[0])
    print("R^2 ajustado:",EstadisticaDelAjuste[1])
    print("chi^2:",EstadisticaDelAjuste[2])
    
    return parameters, paramErrors, EstadisticaDelAjuste

def Ajuste_PA(bin_centers, X, N, errX, p0):
    bc = bin_centers[X>0]
    parameters, covariance = curve_fit(PA, bc, X[bc]/N, sigma = errX[bc]/N, 
                                       p0 = p0, absolute_sigma = True)
    paramErrors = np.sqrt(np.diag(covariance))
    
    print("-------Ajuste PA------")
    
    parametros=['mu', 'p']
    for xe in range(len(parametros)):
        print(parametros[xe],"=",parameters[xe],"+/-",paramErrors[xe])
    
    mu, p = parameters 
    
    ajuste = PA(bc, mu, p)
    
    EstadisticaDelAjuste = stats(ajuste, bc, X[bc]/N, 2)
    
    print("R^2:",EstadisticaDelAjuste[0])
    print("R^2 ajustado:",EstadisticaDelAjuste[1])
    print("chi^2:",EstadisticaDelAjuste[2])
    
    return parameters, paramErrors, EstadisticaDelAjuste

def graph_hist(name, dt, K, ajuste = 'Poisson', labels = ['Datos', 'Poisson'], color = '#ee09d8', 
               vlines = 'True', errors = 'True', alpha = 1):
    tiempos, picos = np.loadtxt(name + '.txt', delimiter=',', unpack=True, skiprows = 1)
    g = contar_eventos(tiempos, dt, .5)
    # print(np.mean(g))
    # K numero de repeticiones del "histograma"
    N = int(len(g)/K)   # numero de repeticiones del experimento (entradas de cada histograma)
    n = int(max(g))     # fotones incidentes
    
    print(' K=',K, ' N=',N)
    
    X, errX = Hist(n, N, K, g)
    
    bin_centers = np.arange(0,n+1,1)
    
    if ajuste == 'Poisson' or ajuste == 'Both':
        parameters, paramErrors, E = Ajuste_poisson(bin_centers, X, N, errX, np.mean(g))
        # mu = np.mean(g)
        p_k = Poisson(bin_centers[X>0], parameters)
        plt.plot(bin_centers[X>0], p_k, marker='s', markersize=5, label='Poisson', color='tab:orange', zorder=10)
    
    if ajuste == 'PA' or ajuste == 'Both':
        p_est = 2/(sdev(g)**2/np.mean(g)+1)
        u_est = np.mean(g)*p_est
        parameters, paramErrors, E = Ajuste_PA(bin_centers, X, N, errX, [u_est, p_est])
        p_k = PA(bin_centers[X>0], parameters[0], parameters[1])
        plt.plot(bin_centers[X>0], p_k, marker='s', markersize=5, label='Polya-Aeppli', color='tab:red', zorder=10)
    
    plt.stairs(X/N, np.hstack((bin_centers,n+1)) -.5 , fill=True, color=color, edgecolor='black',
               baseline=0 , linewidth=2,label=labels[0], alpha = alpha)
    
    if vlines == 'True':
        for i in range(len(bin_centers)):
            xl = bin_centers[i]-.5
            plt.plot([xl,xl],[0,X[i]/N],color='black',linewidth=2,alpha=0.7)
    
    # plt.plot([min(bin_centers), max(bin_centers)],[0,0], color='black', linewidth=2, alpha=1)
    if errors != 'False':
        plt.errorbar(bin_centers, X/N, yerr=errX/N, fmt='none', ecolor='black', capsize=2, elinewidth=1)
    
    if np.sum(X/N) != 1:
        print('Bad norm.', np.sum(X/N))
    
    return N

#%% Ajustes poisson 207Bi y polya fondo

dt = .05
K = 20 # numero de repeticiones del "histograma" (sobre lo cual se promedia) 

plt.figure()
N = graph_hist('fondo_1_resultados', dt, K = 20, ajuste = 'Both', labels = ['Fondo', 'Polya-Aeppli'],
               color = 'tab:blue', vlines = 'False')
plt.title(r'Fondo oscuro (N={N}, K={K})'.format(N = N, K = K),fontsize=14)
plt.xlabel(r'Pulsos detectados en $\Delta t = {}$ s'.format(dt), fontsize=14)
plt.ylabel('Cuentas (Norm.)', fontsize=14)
plt.tight_layout()
plt.xlim(-1, 17)
plt.xticks(np.arange(0,18,2))
plt.legend(fontsize=14)

# plt.savefig('Ajuste_fondo.png', dpi=300)
