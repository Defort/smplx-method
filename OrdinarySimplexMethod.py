#!/usr/bin/env python
# -*- coding: utf-8 -*-
#this script written for solving linear optimization problems
#with simplex method

#import module
import csv;
import numpy as np;
#verbosity on
SHOW_INFO=True;
#aim function coefficients
c=[10,14,12,0,0,0];
#constraints
B=[180,240,210];
#now we defining input matrix of the arguments
INPUT=[[4,2,1,1,0,0],
       [3,1,3,0,1,0],
       [1,2,5,0,0,1]];
#class for solving simplex method problems

class simplexMethod:
    #classification function
    def preproc(INPUT,c):
        if len(c)!=len(INPUT[0]):
            print(len(INPUT[0]));
            print('System cannot solve this problem. Revision is needed.');
            raise ValueError('Cannot perform calculation.')
        if len(c)==len(INPUT[0]):
            print('All clear. Check the basis.');
    #defining basis for the system
    def basis(INPUT):
        SM=[];
        M=np.matrix(INPUT);
        for kit in list(range(len(INPUT[0]))):
            SM.append(sum(M[:,kit]).tolist()[0][0]);
        SM=np.array(SM);
        #defining number of the basis vectors
        BASIS=np.reshape(np.argwhere(SM==1),len(INPUT)).tolist();
        return BASIS;
    #defining indicies of the resolving element
    def resolPosition(INPUT,B,delta):
        #defining resolving column k
        k=np.argmin(delta);
        #defining resolving row r
        n_row=len(INPUT);
        div_r=[];
        for i in list(range(n_row)):
            if INPUT[i][k]>0:
                div_r.append(abs(B[i]/(INPUT[i][k])));
            else:
                div_r.append(1e7);
        r=np.argmin(div_r);
        return r,k;
    def tabulaeINPUT(INPUT,r,k):
        #shape of the matrix
        n_row=len(INPUT);
        n_col=len(INPUT[0]);
        #recalculating our table
        #resolving element
        a_rk=INPUT[r][k];
        #new table
        INPUT_new=[];
        for i in range(n_row):
            for j in range(n_col):
                a_rj=INPUT[r][j];
                a_ik=INPUT[i][k];
                a_ij=INPUT[i][j];
                #triangle rule
                if i!=r:
                    a_ij_new=a_ij-a_rj*a_ik/a_rk;
                elif i==r:
                    a_ij_new=a_ij/a_rk;
                INPUT_new.append(a_ij_new);
        INPUT_new=np.reshape(np.array(INPUT_new),(n_row,n_col));
        return INPUT_new
    def tabulaeDelta(INPUT,delta,r,k):
        #shape of the matrix
        n_row=len(INPUT);
        n_col=len(INPUT[0]);
        #resolving element
        a_rk=INPUT[r][k];
        #new delta
        delta_new=[];
        for i in range(n_col):
            #triangle rule
            a_ri=INPUT[r][i];
            delta_i=delta[i]-delta[k]*a_ri/a_rk;
            delta_new.append(delta_i);
        return delta_new
    def checkDelta(delta):
        pr=[]
        for i in range(len(delta)):
            if delta[i]>=0:
                pr.append(1);
            elif delta[i]<0:
                pr.append(0);
        if sum(pr)==len(delta):
            return True;
        else:
            return False;
    def tabulaeB(INPUT,B,r,k):
        #shape of the matrix
        n_row=len(INPUT);
        n_col=len(INPUT[0]);
        #resolving element
        a_rk=INPUT[r][k];
        #new B
        B_new=[];
        for i in range(n_row):
            #triangle rule
            a_ik=INPUT[i][k];
            if i!=r:
                B_i=B[i]-B[r]*a_ik/a_rk;
            elif i==r:
                B_i=B[i]/a_rk;
            B_new.append(B_i);
        return B_new;
    def tabulaeBasis(BASIS,r,k):
        #assigning k column instead of the r row in BASIS
        i=BASIS[r];
        BASIS[r]=k;
        return BASIS;
    def tabulaeF(INPUT,B,delta,F,r,k):
        #shape of the matrix
        n_row=len(INPUT);
        n_col=len(INPUT[0]);
        #resolving element
        a_rk=INPUT[r][k];
        #new F
        F_new=F-B[r]*delta[k]/a_rk;
        return F_new;
    def solution(INPUT,B,c):
        simplexMethod.preproc(INPUT,c);
        SOLWRITER=csv.writer(open('SOLUTION.csv','a'));
        #first iteration
        BASIS=simplexMethod.basis(INPUT);
        F=0;
        delta=np.array(c)*(-1);
        dchk=simplexMethod.checkDelta(delta);
        SOLWRITER.writerow(INPUT);
        SOLWRITER.writerow(['---']);
        SOLWRITER.writerow(B);
        SOLWRITER.writerow(['---']);
        SOLWRITER.writerow(delta);
        SOLWRITER.writerow(['---']);
        SOLWRITER.writerow(BASIS);
        SOLWRITER.writerow(['---']);
        SOLWRITER.writerow([F]);
        while dchk==False:
            SOLWRITER.writerow(['---']);
            (r,k)=simplexMethod.resolPosition(INPUT,B,delta);
            F=simplexMethod.tabulaeF(INPUT,B,delta,F,r,k);
            B=simplexMethod.tabulaeB(INPUT,B,r,k);
            delta=simplexMethod.tabulaeDelta(INPUT,delta,r,k);
            INPUT=simplexMethod.tabulaeINPUT(INPUT,r,k);
            BASIS=simplexMethod.tabulaeBasis(BASIS,r,k);
            dchk=simplexMethod.checkDelta(delta);
            SOLWRITER.writerow(['r',r,'k',k]);
            SOLWRITER.writerow(INPUT);
            SOLWRITER.writerow(['---']);
            SOLWRITER.writerow(['B=']);
            SOLWRITER.writerow(B);
            SOLWRITER.writerow(['---']);
            SOLWRITER.writerow(['delta=']);
            SOLWRITER.writerow(delta);
            SOLWRITER.writerow(['---']);
            SOLWRITER.writerow(['BASIS=']);
            SOLWRITER.writerow(BASIS);
            SOLWRITER.writerow(['---']);
            SOLWRITER.writerow(['F=']);
            SOLWRITER.writerow([F]);

simplexMethod.solution(INPUT,B,c)
