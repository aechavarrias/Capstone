# -*- coding: utf-8 -*-
"""
Created on Mon Sep 06 21:15:44 2021

Nuevamente modificado por el grupo 5 del ramo Taller de Investigacion Operativa (Capstone).
Se adapto a las especificaciones de la investigacion. 
Se agrego porcentaje de merma, multiples casas y se eliminaron partes sobrantes.
"""
"""
Created on Fri Sep 15 14:23:25 2017

El archivo original fue descargado de http://www.dcc.fc.up.pt/~jpp/mpa/cutstock.py 
y el autor está identificado más abajo
Modificada el 2017-09-15 por J. Vera eliminando la parte heurística y el Bin-Packing para
dejarla como un problema puro de Generación de Columnas
"""
"""
cutstock.py:  use gurobi for solving the cutting stock problem.

The instance of the cutting stock problem is represented by the two
lists of m items of size and quantity s=(s_i) and q=(q_i).

The roll size is B.

Given packing patterns t_1, ...,t_k,...t_K where t_k is a vector of
the numbers of items cut from a roll, the problem is reduced to the
following LP:
    
    minimize   sum_{k} x_k
    subject to sum_{k} t_k(i) x_k >= q_i    for all i
	       x_k >=0			    for all k.

We apply a column generation approch (Gilmore-Gomory approach) in
which we generate cutting patterns by solving a knapsack sub-problem.

Copyright (c) by Joao Pedro PEDROSO and Mikio KUBO, 2010
"""

#from gurobipy import *
from typing import ClassVar
import gurobipy as gp
import numpy as np
from gurobipy import GRB
import random
import math

LOG = True
EPS = 1.e-6



def solveCuttingStock(s, B):
    """solveCuttingStock.
    
    Parameters:
        s - list with item widths
        B - bin capacity

    Returns a solution: list of lists, each of which with the cuts of a roll.
    """
    # print(len(s))
    w = []   # list of different widths (sizes) of items
    q = []   # quantitiy of orders  
    for item in sorted(s):
        if w == [] or item != w[-1]:
            w.append(item)
            q.append(1)
        else:
            q[-1] += 1

    t = []	# patterns
    m = len(w)
    # generate initial patterns with one size for each item width
    for i,width in enumerate(w):
        pat = [0]*m  # vector of number of orders to be packed into one roll (bin)
        pat[i] = int(B/width)
        t.append(pat)

    #if LOG:
        # print ("sizes of orders=", w)
        # print ("quantities of orders=", q)
        # print ("roll size=", B)
        # print ("initial patterns", t)

    iter = 0
    K = len(t)
    master = gp.Model("LP") # master LP problem
    x = {}
    for k in range(K):
        x[k] = master.addVar(obj=1, vtype="I", name="x[%d]"%k)
    master.update()

    orders={}
    for i in range(m):
        coef = [t[k][i] for k in range(K) if t[k][i] > 0]
        var = [x[k] for k in range(K) if t[k][i] > 0]
        orders[i] = master.addConstr(gp.LinExpr(coef,var), ">", q[i], name="Order[%d]"%i)

    master.update()   # must update before calling relax()
    master.Params.OutputFlag = 0 # silent mode
    # master.write("MP" + str(iter) + ".lp")

    while 1:
        iter += 1
        # print("Iteración: ", iter)
        relax = master.relax()
        relax.optimize()
        relax.Params.OutputFlag = 0

        # print ("objective of Master problem:", relax.ObjVal)
        pi = [c.Pi for c in relax.getConstrs()] # keep dual variables

        knapsack = gp.Model("KP")   # knapsack sub-problem
        knapsack.ModelSense=-1   # maximize
        y = {}
        for i in range(m):
            y[i] = knapsack.addVar(obj=pi[i], ub=q[i], vtype="I", name="y[%d]"%i)
        knapsack.update()

        L = gp.LinExpr(w, [y[i] for i in range(m)])
        knapsack.addConstr(L, "<", B, name="width")
        knapsack.update()
        # knapsack.write("KP"+str(iter)+".lp")
        knapsack.Params.OutputFlag = 0 # silent mode
        knapsack.optimize()
        #if LOG:
            # print ("objective of knapsack problem:", knapsack.ObjVal)
        if knapsack.ObjVal < 1+EPS: # break if no more columns
            break

        pat = [int(y[i].X+0.5) for i in y]	# new pattern
        t.append(pat)
        #if LOG:
            # print ("shadow prices and new pattern:")
            #for i,d in enumerate(pi):
            #     print ("\t%5d%12g%7d" % (i,d,pat[i]))
            # print
        
        # add new column to the master problem
        col = gp.Column()
        for i in range(m):
            if t[K][i] > 0:
                col.addTerms(t[K][i], orders[i])
        x[K] = master.addVar(obj=1, vtype="I", name="x[%d]"%K, column=col)
        master.update()   # must update before calling relax()
        # master.write("MP" + str(iter) + ".lp")
        K += 1

    # print ("Fin de la generación de columnas")
    for v in relax.getVars():
            print(v)
    # Finally, solve the IP
    if LOG:
        master.Params.OutputFlag = 1 # verbose mode
    master.optimize()

    if LOG:
        # print 
        # print ("final solution (integer master problem):  objective =", master.ObjVal)
        # print ("Patrones:")
        total_merma = 0
        lista_merma = []
        cortes = 0
        for k in x:
            if x[k].X > EPS:
                print(x[k].X)
                total_escuadria = 0
                lista_cortes = [w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i]) ]
                for i in range(len(lista_cortes)):
                    total_escuadria += lista_cortes[i]
                # print ("Patrón", k,)
                # print ("Cortes:",)
                # print ([w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i]) ],)
                for i in range(m):
                    if t[k][i]>0:
                        cortes += 1
                # print ("--> %d Escuadrías" % int(x[k].X+.5))
                total_merma += (int(x[k].X+.5))*(B-total_escuadria)
                lista_merma.append(((int(x[k].X+.5))*(B-total_escuadria)))
        # print(f"El total de merma generada es: {total_merma} centimetros")

    rolls = []
    for k in x:
        for j in range(int(x[k].X + .5)):
            rolls.append(sorted([w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i])]))
    
    rolls.sort()
    porcentaje_perdida = (total_merma/(B*len(rolls)))*100

    # print(lista_merma)
    # print(len(lista_merma))
    # print("estoy aca")

    return rolls, porcentaje_perdida, lista_merma, cortes, t



def CuttingStockExample1(casas):
    """CuttingStockExample1: create toy instance for the cutting stock problem."""
    # casas = 1      # número de casas
    # w = [113.8, 5.08, 122.8, 22.35, 113.08, 53.56, 273.76, 305.12, 235, 65, 62, 131, 226.63, 152, 251, 246.5, 76.15, 56.92, 213.4, 135.2, 140.7, 48.12, 18, 259.61, 264.02, 5, 1, 36.5, 307.52, 108, 232.38, 221.1, 124.25, 47.69]
    # q = [49, 0, 0, 86, 18, 8, 51, 85, 2, 4, 1, 1, 124, 3, 15, 12, 8, 1, 1, 1, 3, 2, 2, 3, 5, 11, 0, 7, 7, 0, 0, 0, 8, 4] # quantitiy of orders
    
    # w = [113.8, 22.35, 113.08, 53.56, 273.76, 305.12, 235, 65, 62, 131, 226.63, 152, 251, 246.5, 76.15, 56.92, 213.4, 135.2, 140.7, 48.12, 18, 259.61, 264.02,  5, 36.5, 307.52, 124.25, 47.69]
    w = [114, 22, 113, 54, 274, 305, 235, 65, 62, 131, 226, 152, 251, 247, 76, 57, 213, 135, 141, 48, 18, 260, 264,  5, 37, 308, 124, 48]
    q = [49   ,    86,     18,     8,     51,     85,   2,  4,  1,   1,    124,   3,  15,    12,     8,     1,     1,     1,     3,     2,  2,      3,      5, 11,    7,      7,      8,     4] # quantitiy of orders
  
    dict_ = {}
    for elemento in range(len(w)):
        dict_[w[elemento]] = q[elemento]
    for i in range(len(q)):
        q[i] = q[i]*casas
    s = []
    for j in range(len(w)):
        for i in range(q[j]):
            s.append(w[j])

    return s, dict_, q

w = [114, 22, 113, 54, 274, 305, 235, 65, 62, 131, 226, 152, 251, 247, 76, 57, 213, 135, 141, 48, 18, 260, 264,  5, 37, 308, 124, 48]
# w = [113.8, 22.35, 113.08, 53.56, 273.76, 305.12, 235, 65, 62, 131, 226.63, 152, 251, 246.5, 76.15, 56.92, 213.4, 135.2, 140.7, 48.12, 18, 259.61, 264.02,  5, 36.5, 307.52, 124.25, 47.69]


def cuenta_casas(casas_totales, cortes_listos, dict_, patrones_en_orden):
    #la siguiente lista, incluye cuantas casas llevo segun cada corte (el minimo es la cantidad de casas verdaderas)
    lista_de_limitaciones = []
    for i in w:
        if i in cortes_listos.keys():
            limitacion = cortes_listos[i] // dict_[i]
            lista_de_limitaciones.append(limitacion)
        else: # significa que todavia no nos alcanza para hacer una 
            return i 

    a = min(lista_de_limitaciones)
    # for i in range(casas_totales):
    #     if a == i + 1 and f"llevo_{i + 1}_casa" not in patrones_en_orden:
    #         patrones_en_orden.append(f"llevo_{i + 1}_casa")

    # for i in range(casas_totales):
    #     if a == i + 1 and f"llevo_{i + 1}_casa" not in patrones_en_orden:
    #         patrones_en_orden.append(f"llevo_{i + 1}_casa")        
    
    # print(f"EL MINIMO ES: {a}")
    if casas_totales == a: # significa que todos los cortes alcanzan para hacer todas las casas
        # print("terminamos")
        senal_fin = -1
        return senal_fin
    for j in cortes_listos.keys():
        # if "m" in cortes_listos.keys():
            # print("AUXILIOOOOOO cortes_listos")
        # if "m" in dict_.keys():
            # print("AUXILIOOOOOO dict_")
        # print("cagamos")# if j != "m":                                       # OJO GIGANTE, HAY ALGO MALO
        if j in dict_.keys():
            if cortes_listos[j] // dict_[j] == a:
                return j

def elegir_patron(rolls, patrones_en_orden, cortes_listos, corte_limitante):
    if not cortes_listos:
        patron = rolls.pop(0)
        patrones_en_orden.append(patron)
        # print("agregue primer patron")
        for i in patron:
            # cortes_listos.append(i)
            if i in cortes_listos.keys():
                cortes_listos[i] += 1
            else:
                cortes_listos[i] = 1
    else:
        for i in range(len(rolls)):
            poped = False
            if corte_limitante in rolls[i]:
                patron = rolls.pop(i)
                poped = True
                patrones_en_orden.append(patron)
                # print(f"agregue patron nuevo {patron}")
                for k in patron:
                    # cortes_listos.append(k)
                    if k in cortes_listos.keys():
                        cortes_listos[k] += 1
                    else:
                        cortes_listos[k] = 1
            if poped:
                break

def ordenar_patrones(lista_desordenada, casas, cortes_listos, dict_):
    lista_ordenada = []
    for patron in range(len(lista_desordenada)):
        # print(type(cortes_listos))
        corte_limitante = cuenta_casas(casas, cortes_listos, dict_, lista_ordenada)
        if corte_limitante == -1:
            # print("estoy cortando")
            break
        else:
            elegir_patron(lista_desordenada, lista_ordenada, cortes_listos, corte_limitante)
    return lista_ordenada

def mejor_largo(largo_maximo):
    mermas = []
    for i in range(largo_maximo):
        rolls, merma = solveCuttingStock(s, i + 308)
        mermas.append(merma)
    a = min(mermas)
    b = mermas.index(a)
    c = b + 308
    return a, c

def patrones_por_casa(lista_ordenada):
    contador = 0
    lista_patrones_por_casa = []
    # lista_ordenada.append("llevo_50_casa")
    for i in lista_ordenada:
        tipo = isinstance(i, str)
        if tipo == False:
            contador += 1
        else:
            lista_patrones_por_casa.append(contador)
            contador = 0
    return lista_patrones_por_casa

def ordenar_mermas(patrones, largo):
    mermas = []
    for tabla in patrones:
        usado = sum(tabla)
        merma = largo - usado
        mermas.append(merma)
    mermas.sort()
    return mermas

def funcion(lista_ordenada, lista_mermas, largo, casas, dict_, escuadrias_utilizadas, lista_completos):
    # lista_completos = {}
    lista_faltantes = []
    lista_incompletos = []
    trozo_malo = 0
    largo_escuadria = 0
    for escuadria in lista_ordenada:
        for trozo in escuadria:
            tipo = isinstance(trozo, str)
            if tipo == False:
                # print(trozo)
                lista_faltantes.append(float(trozo))
    # print(f"LA CANTIDAD TOTAL DE TROZOS ES: {len(lista_faltantes)}\n")
    for escuadria in lista_ordenada:

        # if escuadria[-2] == "m":
        # print(escuadria)
        # if sum(escuadria) == 320:
        #     print("CHECKPOINT 2")
        #     escuadrias_utilizadas += 1
        #     largo_escuadria = largo
        # else:
        #     print("CHECKPOINT 1")
        #     largo_escuadria = sum(escuadria[:-1])
        #     lista_mermas.append(escuadria[-1])
        #     escuadria = escuadria[:-2]
        # for trozo in escuadria:
        # if trozo == "m":
        #     print("CHECKPOINT 1")
        #     largo_escuadria = sum(escuadria[:-2])
        #     lista_mermas.append(escuadria[-1])
        #     escuadria = escuadria[:-2]
        #     break
        # if largo_escuadria == 0:
        #     print("CHECKPOINT 2")
        #     escuadrias_utilizadas += 1
        print(f"la suma de la escuadria es: {sum(escuadria)}")
        print(f"la escuadria es {escuadria}")
        print(f"el largo es {largo}")
        if sum(escuadria) <= largo:
            largo_escuadria = largo
        else:
            print("entre al else")
            largo_escuadria = escuadria[-1] / largo
            print(f"la merma es {merma}")
            escuadria = escuadria[:-1]

        for trozo in escuadria:
            indicador = False
            if trozo > largo_escuadria:  #and (trozo - largo_escuadria) > 1
                # print(trozo, ">", largo_escuadria)
                lista_incompletos.append(trozo)
                lista_mermas.append(largo_escuadria)
                # print("Me quede sin espacio en la escuadria")
                continue
            probabilidad = random.random()
            # print(f"la probabilidad es {probabilidad}")
            if 0.25 < probabilidad < 0.75:# prob salga bien cortado
                lista_faltantes.remove(trozo)
                # lista_completos.append(trozo)
                if trozo in lista_completos.keys():
                    lista_completos[trozo] += 1
                else:
                    lista_completos[trozo] = 1
                largo_escuadria -= trozo
            elif 0 < probabilidad < 0.05: # -2
                trozo_malo = trozo - 2
                lista_mermas.append(trozo_malo)
                lista_mermas.append(largo_escuadria) # manejar
                lista_mermas.sort()
                
                rolls, lista_mermas, lista_faltantes = optimizar(lista_faltantes, largo, lista_mermas)
                lista_ordenada = ordenar_patrones(rolls, casas, lista_completos, dict_)
                lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas = funcion(lista_ordenada, lista_mermas, largo, casas, dict_, escuadrias_utilizadas, lista_completos)
                return lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas
                # for merma in lista_mermas:
                #     if merma >= trozo and indicador == False:
                        # merma = merma - trozo
                        # lista_faltantes.remove(trozo)
                        # lista_completos.append(trozo)
                        # lista_mermas.append(trozo_malo)
                        # lista_mermas.sort()
                        # largo_escuadria -= (trozo_malo)
                        # indicador = True
                # if indicador == False:
                #     lista_incompletos.append(trozo)
                #     lista_mermas.append(trozo_malo)
                #     print("no cabia en ninguna otra merma (-2 cm)")
            elif 0.05 < probabilidad < 0.25: # -1
                trozo_malo = trozo - 1
                lista_mermas.append(trozo_malo)
                lista_mermas.append(largo_escuadria) # manejar
                lista_mermas.sort()
                
                rolls, lista_mermas, lista_faltantes = optimizar(lista_faltantes, largo, lista_mermas)
                # print(f"AYUDAAAAAAAAAAAAA {lista_ordenada} \n {lista_completos}")
                # if "m" in lista_completos:
                #     print("estamos cagados")
                lista_ordenada = ordenar_patrones(rolls, casas, lista_completos, dict_)
                # print(f"AYUDAAAAAAAAAAAAA {lista_ordenada} \n {lista_completos}")
                lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas = funcion(lista_ordenada, lista_mermas, largo, casas, dict_, escuadrias_utilizadas, lista_completos)
                return lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas
                # for merma in lista_mermas:
                #     if merma >= trozo and indicador == False:
                #         merma = merma - trozo
                #         lista_faltantes.remove(trozo)
                #         lista_completos.append(trozo)
                #         lista_mermas.append(trozo_malo)
                #         lista_mermas.sort()
                #         largo_escuadria -= (trozo_malo)
                #         indicador = True
                # if indicador == False:
                #     lista_incompletos.append(trozo)
                #     lista_mermas.append(trozo_malo)
                #     print("no cabia en ninguna otra merma (-1 cm)")
            elif 0.75 < probabilidad < 0.95: # 1
                trozo_malo = trozo + 1
                largo_escuadria -= (trozo_malo)
                lista_faltantes.remove(trozo)
                # lista_completos.append(trozo)
                if trozo in lista_completos.keys():
                    lista_completos[trozo] += 1
                else:
                    lista_completos[trozo] = 1
                lista_mermas.append(1)
                lista_mermas.sort()
            elif 0.95 < probabilidad < 1:   # 2
                trozo_malo = trozo + 2
                largo_escuadria -= (trozo_malo)
                lista_faltantes.remove(trozo)
                # lista_completos.append(trozo)
                if trozo in lista_completos.keys():
                    lista_completos[trozo] += 1
                else:
                    lista_completos[trozo] = 1
                lista_mermas.append(2)
                lista_mermas.sort()
    # print(f"EN LA LISTA FALTANTES HAY {len(lista_faltantes)} TROZOS")
    return lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas

def optimizar(restantes, largo, lista_mermas):
    new_rolls = []
    new_lista_merma = []
    new_restantes = []
    for merm in lista_mermas:
        lista = []
        largo_patron = 0
        # if merm >= 308:
        #     rolls, merma, lista_merma, cortes, t = solveCuttingStock(restantes, merm)
        #     patron_elegido = elegir_mejor(rolls)
        #     for corte in patron_elegido:
        #         largo_patron += corte
        #         restantes.remove(corte)
        #     # new_lista_merma.append(merm-largo_patron)
        #     # patron_elegido.append("m")
        #     patron_elegido.append(merm-largo_patron)
        #     new_rolls.append(patron_elegido)
        if 100 < merm < 308: # ojo
            lista = cortes_menores_a(restantes, merm)
            rolls, merma, lista_merma, cortes, t = solveCuttingStock(lista, merm)
            patron_elegido = elegir_mejor(rolls)
            for corte in patron_elegido:
                restantes.remove(corte)
                largo_patron += corte
            new_lista_merma.append(merm-largo_patron)
            # patron_elegido.append("m")
            # patron_elegido.append(merm-largo_patron)
            extra = merm * largo
            patron_elegido.append(extra)
            new_rolls.append(patron_elegido)
        else:
            new_lista_merma.append(merm)
    rolls, merma, lista_merma, cortes, t = solveCuttingStock(restantes, largo)
    # print("optimizo 2.5")
    for patron in rolls:
        # print(patron)
        new_rolls.append(patron)
    # for excedente in lista_merma:
    #     new_lista_merma.append(excedente)
    for patron in new_rolls:
        largo_patron = sum(patron)
        merma_patron = largo - largo_patron
        new_lista_merma.append(merma_patron)
        for cut in patron:
            new_restantes.append(cut)
    return new_rolls, new_lista_merma, new_restantes
    # return

def cortes_menores_a(lista_mermas, largo):
    lista_return = []
    for merma in lista_mermas:
        if merma < largo:
            lista_return.append(merma)
    return lista_return

def elegir_mejor(rolls): # Cambiar  (elegir el con menos merma) (elegir el con menos cortes o cortes mas largos) (elegir el con menos cortes)
    # print(rolls)
    patron_elegido = [0]
    for patron in rolls:
        if len(patron) > len(patron_elegido):
            patron_elegido = patron
    return patron_elegido

if __name__ == "__main__":

    # escuadrias_utilizadas = 0

    # cortes_listos = {}
    # lista_mermas = []
    lista_escuadrias = []
    # casas = 1
    # largo = 320
    for i in range(1):
        escuadrias_utilizadas = 0
        cortes_listos = {}
        lista_mermas = []
        # lista_escuadrias = []
        casas = 1
        largo = 320
        s, dict_, q = CuttingStockExample1(casas)
        rolls, merma, lista_merma, cortes, t = solveCuttingStock(s, largo)
        print(rolls)
        # for escuadria in rolls:
        #     largo_sin_merma = sum(escuadria)
        #     merma_escuadria = largo - largo_sin_merma
        #     escuadria.append(merma_escuadria)
        # rolls2, lista_merma2 = optimizar(s, largo, lista_mermas = [])
        '''
        # rolls, porcentaje_perdida, lista_merma, cortes, t

        # largo_maximo = 800
        # porcentaje_merma, largo_elegido = mejor_largo(largo_maximo)

        # print (len(rolls), "Escuadrías:")
        '''
        # print(f"los rolls son los siguientes: {rolls}")
        lista_mermas = ordenar_mermas(rolls, largo)
        # print(f"la lista es la siguiente: {lista_mermas}")
        lista_ordenada = ordenar_patrones(rolls, casas, cortes_listos, dict_)

        '''
        lista_patrones_por_casa = patrones_por_casa(lista_ordenada)
        '''

        # print(type(lista_ordenada[1][1]))
        lista_completos = {}
        lista_completos, lista_incompletos, escuadrias_utilizadas, lista_mermas = funcion(lista_ordenada, lista_mermas, largo, casas, dict_, escuadrias_utilizadas, lista_completos)

        # print(f"en total hay {escuadrias_utilizadas} escuadrias")
        # print(f" la lista completos: {lista_completos} \n y la lista incompletos: {lista_incompletos} \n y la lista mermas {lista_mermas}")

        total = 0

        for key in lista_completos:
            total += key * lista_completos[key]
        
        for merma in lista_mermas:
            total += merma
        escuadrias = total/320
        lista_escuadrias.append(escuadrias)
    print(lista_escuadrias)
    # promedio = sum(lista_escuadrias)/50
    # print(promedio)
    # print(total)
    # print(total/320)

    # print(dict_)
    # print(f"el largo de la lista de trozos completos es {len(lista_completos)}, el de la lista de trozos incompletos es {len(lista_incompletos)} y su suma es {len(lista_completos) + len(lista_incompletos)}")

    # print(lista_incompletos)

    # rolls, merma, lista_merma, cortes, t = solveCuttingStock(lista_incompletos, largo)

    # print(rolls)
    # print(len(rolls))

 
"""
pescamos el patron, vemos si hay error (20% error por patron, que debiese ser por corte)
en caso de error (60% no se equivoque, 30% +/- 1, 10% +/-2), el corte con error mas grande, se corta bien (supuesto de que siempre se corta bien)
el corte con error mas chico se intenta agregar a algun otro patron con merma disponible, en caso de no lograrlo,
se agrega a los cortes faltantes para las casas y se reoptimiza el problema.

"""




    

    
    