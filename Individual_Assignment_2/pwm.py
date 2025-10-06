#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 11:16:02 2024

@author: inestavares
"""
import numpy as np
from Bio.Seq import Seq
from Bio import motifs
class PWM:
    def __init__ (self, lst_of_words=None, bio_type='DNA',alphabet='ACTG',
                  pwm=None,pseudo=0.01): 
        self.lst_of_words = lst_of_words if lst_of_words else []
        self.bio_type = bio_type
        self.alphabet = alphabet
        self.pwm_size = len(alphabet)
        self.pseudo=pseudo
        if pwm:
            self.pwm = pwm
        else:
            pwm_size = 4
            self.pwm = [[0] * pwm_size for _ in range(pwm_size)]
        

    def arguments(self):
        print('List of Motif Words: ')
        if self.lst_of_words!= None:
            l=len(self.lst_of_words)
            for i in range(l):
                print(self.lst_of_words[i])
        else:
            print('')
        print('\n')
        print('Type of Sequence:',self.bio_type,'\n')
        print('Alphabet:', self.alphabet,'\n')
        print('Frequency Matrix:')
        self.r=len(self.pwm)
        for i in range(self.r):
            print(self.alphabet[i],':',self.pwm[i])



    def add_pseudocounts(self):
        self.r=len(self.pwm)
        self.c=len(self.pwm[0])
        self.pwm_pseudo = [[elem + self.pseudo for elem in row] for row in self.pwm]
        return(self.pwm_pseudo)
    
            
    def print_pwm(self):
        self.pwm_pseudo=self.add_pseudocounts()
        self.r=len(self.pwm_pseudo)
        self.c=len(self.pwm_pseudo[0])
        print('Matrix PWM with a pseudo count of ',self.pseudo)
        for i in range(self.r):
            print(self.alphabet[i],':',self.pwm_pseudo[i])

    
    def informationContent(self):
        self.pwm_pseudo=self.add_pseudocounts()
        self.r=len(self.pwm_pseudo)
        self.c=len(self.pwm_pseudo[0])
        I=[]
        for p in range(self.c):
            Ii=2
            for a in range(self.r):
                prob = self.pwm_pseudo[a][p]
                Ii += prob * np.log(prob) / np.log(2)
            I.append(Ii)
        self.I=I
        print('Information Content of the Motif.')
        for a in range(self.r):
            print(self.alphabet[a],':',self.I[a])

    

def test_default():
    print('Test with the default values.')
    pwm_default = PWM()
    print('\n')
    pwm_default.arguments()
    print('\n')
    pwm_default.print_pwm()
    print('\n')
    pwm_default.informationContent()



def test():
    m1='TATAAAA'
    m2='TATAAAT'
    m3='TATATAT'
    m4='TATAAGG'
    m5='TATAATG'
    m6='CATAAAA'
    m7='CCTATAA'
    m8='TATAATC'
    l=[m1,m2,m3,m4,m5,m6,m7,m8]
    b_type='DNA'
    alph='ACTG'
    p=[[0,0.875,0,1,0.75,0.625,0.375],
         [0.25,0.125,0,0,0,0,0.25],
         [0.75,0,1,0,0.25,0.25,0.25],
         [0,0,0,0,0,0.25,0.125]
         ]
    print('Test with the non default values.')
    print('\n')
    a=PWM(l, b_type, alph, p)
    a.arguments()
    print('\n')
    a.print_pwm()
    print('\n')
    a.informationContent()
    
    


        
if __name__ == "__main__":
    print('\n')
    test_default()
    print('\n')
    print('\n')
    test()
    
    
            
        
            
            
        
    
    




