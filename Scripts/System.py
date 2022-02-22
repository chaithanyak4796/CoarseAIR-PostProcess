#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:15:49 2022

@author: chaithanya
"""
import numpy as np
import sys

au_s = 2.4188843265E-17;
Bo_angs = 0.529177

m_O = 29148.94559
m_N = 25526.04298

class System:
    def __init__(self, sys_name):
        
        self.sys_name = sys_name
        
        if(sys_name == 'O3'):
            self.m1 = m_O
            self.m2 = m_O
            self.m3 = m_O
            
            self.num_reac   = 1
            self.arr_matrix = [[16.0, 17.0, 19.0, 32.0, 33.0, 35.0, 48.0, 49.0, 51.0]]
            
        elif(sys_name == 'N3'):
            self.m1 = m_N
            self.m2 = m_N
            self.m3 = m_N
            
            self.num_reac   = 1
            self.arr_matrix = [[16.0, 17.0, 19.0, 32.0, 33.0, 35.0, 48.0, 49.0, 51.0]]
            
        elif(sys_name == 'N2O'):
            self.m1 = m_N
            self.m2 = m_N
            self.m3 = m_O
            
            self.num_reac   = 2
            self.arr_matrix = [[16.0, 17.0, 19.0], [32.0, 33.0, 35.0, 48.0, 49.0, 51.0]]
            
        else:
            print(" ERROR : Invalid sys_name. Exiting.")
            sys.exit()
            