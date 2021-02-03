# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 19:34:30 2021

@author: caioh
"""

import numpy as np
from correlacoes_empiricas3 import f_empiricas

####inputs###########
########### Do reator ####################
P     = 7       #MPa   (pressão de reação)
#D     = .78     # cm  (Diâmetro do reator)
#h     =         # cm (comprimento do reator) 
#Area  = .478    # cm² (Area da secao) 
T     = 613.15  #K     (Temperatura)


rho0  = .89  #g/cm³ (densidade do óleo)
API   = 30      #°API
Vgas  = 10      # a.u   (razao H2/carga)
F     = .00235  # cm³/s (vazão de carga)

########### Do catalisador ####################
dp    = .035 # cm (diâmetro da partícula) #dissertação ana 2012
#dpe   = .02 # cm (diametro da particula equivalente do catalisador mais inerte)
E     = .53  #    (porosidade do catalisador)
e     = .52  #    (fracao de vazios no leito)
qci   = .55  #    (razao do leito catalitico diluido por inertes)

rhoB  = .75  # g/cm³ (densidade do bulk catalitico)
rhos  = .75  # g/cm³ (densidade do catalisador)
nHDT  = .7   # (fator efetividade do catalisador para um processo de HDT)


eps_B = .4   # a.u (porosidade do leito)



########Cálculo dos parâmetros empíricos#################
resultado = np.zeros(11)
resultados  = f_empiricas(P,T, rho0, API, Vgas, F, dp, eps_B)


uG = resultados[0]
uL = resultados[1] 

k_L_H2al  = resultados[2]
k_L_H2Sal = resultados[3]
k_L_NH3al = resultados[4]
 
k_s_H2  = resultados[5]
k_s_H2S = resultados[6]
k_s_NH3 = resultados[7]

a_s = resultados[8]

k_s_s = resultados[9]
k_s_N = resultados[10]
#########################################################

