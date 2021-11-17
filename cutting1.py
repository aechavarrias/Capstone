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
import gurobipy as gp
import numpy as np
from gurobipy import GRB

LOG = True
EPS = 1.e-6

def solveCuttingStock(s, B):
    """solveCuttingStock.
    
    Parameters:
        s - list with item widths
        B - bin capacity

    Returns a solution: list of lists, each of which with the cuts of a roll.
    """
    print(len(s))
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

    if LOG:
        print ("sizes of orders=", w)
        print ("quantities of orders=", q)
        print ("roll size=", B)
        print ("initial patterns", t)

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
        print("Iteración: ", iter)
        relax = master.relax()
        relax.optimize()

        print ("objective of Master problem:", relax.ObjVal)
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
        if LOG:
            print ("objective of knapsack problem:", knapsack.ObjVal)
        if knapsack.ObjVal < 1+EPS: # break if no more columns
            break

        pat = [int(y[i].X+0.5) for i in y]	# new pattern
        t.append(pat)
        if LOG:
            print ("shadow prices and new pattern:")
            for i,d in enumerate(pi):
                print ("\t%5d%12g%7d" % (i,d,pat[i]))
            print
        
        # add new column to the master problem
        col = gp.Column()
        for i in range(m):
            if t[K][i] > 0:
                col.addTerms(t[K][i], orders[i])
        x[K] = master.addVar(obj=1, vtype="I", name="x[%d]"%K, column=col)
        master.update()   # must update before calling relax()
        # master.write("MP" + str(iter) + ".lp")
        K += 1

    print ("Fin de la generación de columnas")
    for v in relax.getVars():
            print(v)
    # Finally, solve the IP
    if LOG:
        master.Params.OutputFlag = 1 # verbose mode
    master.optimize()

    if LOG:
        print 
        print ("final solution (integer master problem):  objective =", master.ObjVal)
        print ("Patrones:")
        total_merma = 0
        for k in x:
            if x[k].X > EPS:
                total_escuadria = 0
                lista_cortes = [w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i]) ]
                for i in range(len(lista_cortes)):
                    total_escuadria += lista_cortes[i]
                print ("Patrón", k,)
                print ("Cortes:",) 
                print ([w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i]) ],)
                print ("--> %d Escuadrías" % int(x[k].X+.5))
                total_merma += (int(x[k].X+.5))*(B-total_escuadria)
        print(f"El total de merma generada es: {total_merma} centimetros")

    rolls = []
    for k in x:
        for j in range(int(x[k].X + .5)):
            rolls.append(sorted([w[i] for i in range(m) if t[k][i]>0 for j in range(t[k][i])]))
    
    # print("-------------------------------------------------------------------------------------------")
    # print(rolls)
    
    rolls.sort()
    # print(rolls)
    # print("-------------------------------------------------------------------------------------------")
    porcentaje_perdida = (total_merma/(B*len(rolls)))*100

    return rolls, porcentaje_perdida



def CuttingStockExample1(casas):
    """CuttingStockExample1: create toy instance for the cutting stock problem."""
    B = 320         # roll width (bin size) 
    # casas = 1       # número de casas
    w = [113.8, 5.08, 122.8, 22.35, 113.08, 53.56, 273.76, 305.12, 235, 65, 62, 131, 226.63, 152, 251, 246.5, 76.15, 56.92, 213.4, 135.2, 140.7, 48.12, 18, 259.61, 264.02, 5, 1, 36.5, 307.52, 108, 232.38, 221.1, 124.25, 47.69]
    q = [49, 0, 0, 86, 18, 8, 51, 85, 2, 4, 1, 1, 124, 3, 15, 12, 8, 1, 1, 1, 3, 2, 2, 3, 5, 11, 0, 7, 7, 0, 0, 0, 8, 4] # quantitiy of orders
    for i in range(len(q)):
        q[i] = q[i]*casas
    s = []
    for j in range(len(w)):
        for i in range(q[j]):
            s.append(w[j])

    return s,B

def cuenta_casas(lista):
    pass

def elegir_patron():
    pass

def desrelajar(s, B, escuadrias, new_rolls):
    # escuadrias = 0
    rolls, merma = solveCuttingStock(s,B)
    patron = rolls.pop(0)
    for trozo in patron:
        if trozo in s:
            s.remove(trozo)
            if trozo in new_rolls:
                new_rolls[trozo] += 1
            else:
                new_rolls[trozo] = 1
        else:
            escuadrias, new_rolls = desrelajar(s,B,escuadrias, new_rolls)
    escuadrias += 1
    while s:
        rolls, merma = solveCuttingStock(s,B)
        patron = rolls.pop(0)
        print(f"la lista de restantes es {s} \ny el patron es {patron}")
        for trozo in patron:
            print(trozo)
            if trozo in s:
                s.remove(trozo)
                if trozo in new_rolls:
                    new_rolls[trozo] += 1
                else:
                    new_rolls[trozo] = 1
            else:
                escuadrias, new_rolls = desrelajar(s,B,escuadrias, new_rolls)
        escuadrias += 1
    return escuadrias, new_rolls


if __name__ == "__main__":

    cortes_listos = []
    dic = {}
    casas = 1

    new_rolls = {}
    escuadrias = 0

    s,B = CuttingStockExample1(casas)
    escuadrias, new_rolls = desrelajar(s,B, escuadrias, new_rolls)
    print(f"las escuadrias son: {escuadrias}")

    print(f"el diccionario es: {new_rolls}")
    # print("\n\n\nCutting stock problem:")
    # rolls, merma = solveCuttingStock(s,B)
    # largo_total = rolls * B

    # print(rolls)

    # print(len(rolls), "Escuadrías:")
    # print(f"El porcentaje de merma generada es {merma}%")

    # for escuadria in rolls:
    #     for trozo in escuadria:
    #         if trozo in dic:
    #             dic[trozo] += 1
    #         else:
    #             dic[trozo] = 1

    # print(dic)
    