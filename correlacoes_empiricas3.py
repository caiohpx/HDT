# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 18:35:28 2021

@author: caioh
"""

import numpy as np


####inputs###########
P     = 7       #MPa   (pressão)
T     = 613.15  #K     (Temperatura)
rho0  = .89  #g/cm³ (densidade do óleo)
API   = 30     #°API
#RH2Hc = 827.234 # a.u   (razao H2/carga)
Vgas  =  10        # vazão do gás em NL/h 
F     = .00235  # cm³/s (vazão de carga)


dp    = .035 # cm  (diâmetro da partícula)
eps_B = .4   # a.u (porosidade do leito)



def f_empiricas(P,T, rho0, API, Vgas, F, dp, eps_B):             # imput ==> P em MPa , T em K ,  desnidade do óleo rhoL em g/cm³, API,  vazão do gás em NL/h, cm³/s vazao de carga ; 0,141mL/s, cm diâmetro da partícula, porosidade do leito


    ####      conversões   #####
    T_R        = T*1.8         # °R      Obs: 1K      = 1,8     °R (graus Rakine)
    rho0_lbft3 = rho0*62.428   # lb/ft3  Obs: 1g/cm³  = 62,428  lb/ft3     Obs: a densidade em g/cm³ tem que ser a 20°C
    P_psi      = P*145.038     # psi     Obs: 1 MPa   = 145,038 psi
    
    #print(rho0_lbft3)
    
    
    
    #### Densidade do óleo #####
    Termo1     = .167+(16.181*(10**(-.0425*rho0_lbft3)))
    Termo2     = .299+(263*(10**(-.0603*rho0_lbft3)))
    deltarhoP  = Termo1*(P_psi/1000) - .01*Termo2*((P_psi/1000)**2)
    
    
    
    Termo3     = .0133+(152.4*((rho0_lbft3+deltarhoP)**-2.45))
    Termo4     = 8.1e-6 - .0622*(10**(-.764*(rho0_lbft3+deltarhoP)))
    deltarhoT  = Termo3*(T_R-520) - Termo4*((T_R-520)**2)
    
    rhoL_lbft3 = rho0_lbft3 + deltarhoP - deltarhoT
    
    rhoL       = rhoL_lbft3/62.428
    
    
    #### Solubilidade dos gases #####
    def lambdi(T, rhoL):
        
        lambdH2  = -.559729-((.42947e-3)*(T-273.15)) + (3.07539e-3)*((T-273.15)/rhoL) + (1.94593e-6)*((T-273.15)**2) + (.835783/(rhoL**2)) #  L/(Kg.MPa) Solubilidade de H2. T em K e densidade (rho_L) na temperatura de reação
        lambdH2S = np.exp(3.367-(0.00847*(T-273.15)))                       # L/(Kg.MPa) Solubilidade de H2S
        lambdNH3 = (1/((8.552e-2)+((2.233e-6)*(T**2.79))))*1000    # Chacón et. al 2012
        return lambdH2, lambdH2S, lambdNH3 #L/(Kg.MPa) 
    
    lambdH2, lambdH2S, lambdNH3 = lambdi(T, rhoL)
    
    
    
    
    def Henry(lambdi):
        nuN  = 22.413  # L/mol Volume molar dos gases nas condições padrão
        
        return (1e3*nuN/(rhoL*lambdi))
    
    HH2, HH2S, HNH3 = Henry(lambdH2), Henry(lambdH2S), Henry(lambdNH3)
    
    #############conferido#################
    
    
    
    ### Viscosidade dinâmica do líquido#############
    a    = (10.313*(np.log10((T_R - 460)))) - 36.447        # dinâmica do líquido (a.u) 
    muL = (3.141e10)*((T_R-460)**-3.444)*((np.log10(API))**a) # viscoisdade (cpoise)
    
    
    
    
    
    ########Volume molar no líquido (nu_i)######
    TMeABP = 680.15 # ponto de ebulicao medio em Kelvin";#artigo modeling of a three-phase reactor Chacon et al#
    MML    = .7     # lb/mol   MASSA MOLAR DO DIESEL em lb/mol ref.: artigo modeling of a three-phase reactor Chacon et al / C15H28 + contaminantes#
    nuc    = (MML*(7.5214e-3)*((TMeABP*1.8)**.2896)*(rho0_lbft3**(-.7666))) * (28316.8) # cm³/mol Volume molar crítico. o último termo passa de ft³ para cm³. o 1.8 passa para °R
    
    def VolMol(nuc_i):
            
            
        nu_i  = 0.285*(nuc_i**1.048) # volume molar do líquido #cm³/mol
            
        return nu_i
    
    
    nuL, nu_S, nu_N  = VolMol(nuc), VolMol(nuc), VolMol(nuc) # cm³/mol
    nu_H2, nu_H2S, nu_NH3 = 616, 571, 576
    
    
    ##########DIFUSIVIDADE DO LÍQUIDO#############
        
    def Difusividade(nu_i):   # cm²/s Difudividade molecular de H2  na mistura 
        
        Difus_i = (8.93e-8)*((nuL**.267)/(nu_i**.433))*(T/(muL)) 
        
        return Difus_i 
    
    Difus_H2, Difus_H2S, Difus_NH3, Difus_S, Difus_N = Difusividade(nu_H2), Difusividade(nu_H2S), Difusividade(nu_NH3), Difusividade(nu_S), Difusividade(nu_N)
    
    
    
    
    
    ##########COEFICIENTES DE TRANSFERÊNCIA DE MASSA GÁS-LÍQUIDO###############
    Area  =  .478     # cm²
    e     =  .52       # fracoes de vazios no leito
    #Vgas  =  10        # vazão do gás em NL/h (input)
    # F     = 0.00235  # cm³/s vazao de carga ; 0,141mL/s (input)
    Fg    = .1*Vgas*T/(273.15*.277*P) # 1.15 cm³/s Vazão do gás; 
    
    
    uG = Fg/(Area*e)        # cm/s velocidade superficial do gás
    uL = F/(Area*e)         # cm/s velocidade superficial do liquido
    GL = uL*(rhoL_lbft3*453.592/28316.8)            # g.s-1.cm-2 fluxo de massa da fase liquida
    
    
    
    def coefGL(Difus_i):
        
        k_ial = (Difus_i*7)*(((GL/(muL))**.4))*(abs(((muL*.01)/(rhoL*Difus_i))**.5))  #1/s 0.000191894    # 1/s  o .01 é para deixar o ultimo termo entre parenteses admensional, passando de cp para g/(cm*s)
       
        return k_ial
    
    k_L_H2al, k_L_H2Sal, k_L_NH3al  = coefGL(Difus_H2) , coefGL(Difus_H2S), coefGL(Difus_NH3)
    
    
    
    
    ##########COEFICIENTES DE TRANSFERÊNCIA DE MASSA LÍQUIDO-SÓLIDO###############
    
    dpe = .02               #diametro da particula equivalente do catalisador mais inerte (input)
    a_s = 6*(1-e)/dpe       #1/cm 102.857 1/cm Área superfícial específica (área específica do sólido)
    
    
    def coefLS(Difus_i):            ##########PAREI AQUI
        #rhoL_lbft3 = 46.1072
        Termoa = 1.8*(Difus_i*a_s)
        Termob = abs((GL/(a_s*muL*.01))**.5)
        Termoc = abs((muL*.01/(rhoL*Difus_i))**(1/3))
    
        return Termoa*Termob*Termoc               
    
    k_s_H2, k_s_H2S,  k_s_NH3, k_s_s, k_s_N = coefLS(Difus_H2), coefLS(Difus_H2S), coefLS(Difus_NH3), coefLS(Difus_S), coefLS(Difus_N)

    return np.array([uG, uL, k_L_H2al, k_L_H2Sal, k_L_NH3al, k_s_H2, k_s_H2S,k_s_NH3, a_s, k_s_s, k_s_N])
  

