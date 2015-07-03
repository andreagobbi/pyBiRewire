cimport cython
import numpy as np
cimport numpy as np
from math import *
import igraph as i
import os
import csv

cdef extern from "lib/BiRewire.h":
    double similarity_undirected(unsigned short *m,unsigned short *n,size_t ncol,size_t nrow,size_t e)
    double similarity(unsigned short *m,unsigned short *n,size_t ncol,size_t nrow,size_t e )
    size_t analysis_ex(unsigned short *incidence,size_t ncol, size_t nrow,double *scores,size_t step,size_t max_iter,size_t verbose,size_t MAXITER)
    size_t analysis(unsigned short *incidence,size_t ncol, size_t nrow,double *scores,size_t step,size_t max_iter,size_t verbose)
    size_t rewire_bipartite(unsigned short *matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose)
    size_t rewire_bipartite_ex(unsigned short *matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER)
    size_t rewire_sparse_bipartite_ex(size_t *fro,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose,size_t MAXITER)
    size_t rewire_sparse_bipartite(size_t *fro,size_t *to,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose)

    size_t analysis_undirected(unsigned short *incidence,size_t ncol, size_t nrow,double *scores,size_t step,size_t max_iter,size_t verbose)
    size_t analysis_undirected_ex(unsigned short *incidence,size_t ncol, size_t nrow,double *scores,size_t step,size_t max_iter,size_t verbose,size_t MAXITER)
    size_t rewire(unsigned short *matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose)
    size_t rewire_ex(unsigned short *matrix,size_t ncol, size_t nrow,size_t max_iter,size_t verbose,size_t MAXITER)
    #size_t rewire_sparse_ex(size_t *fro,size_t *to,size_t *degree,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose,size_t MAXITER)
    #size_t rewire_sparse(size_t *fro,size_t *to,size_t *degree,size_t nc,size_t nr,size_t max_iter,size_t ne,size_t verbose)

def c_rewire_sparse_bipartite(np.ndarray left,np.ndarray right,N=-1, verbose=1,  MAXITER=10, accuracy=0.00005,exact=True):

    cdef size_t e,nc,nr,t
    e= len(left)
    nc,nr= len(np.unique(left)),len(np.unique(right))
    t=nc*nr
    if N ==-1:
        N=ceil((e*(1-e/t)) *log((1-e/t)/accuracy) /2  )  
    if exact:
        return N,<int>rewire_sparse_bipartite_ex(<size_t*>np.PyArray_DATA(left),<size_t *>np.PyArray_DATA(right),nr, nc, N, e, verbose, N*MAXITER)
    else:
        return N,<int>rewire_sparse_bipartite(<size_t*>np.PyArray_DATA(left),<size_t *>np.PyArray_DATA(right),nr, nc, N, e, verbose)


def c_rewire_bipartite(np.ndarray incidence,N=-1, verbose=1,  MAXITER=10, accuracy=0.00005,exact=True):
    cdef size_t e,nc,nr,t
    nc,nr= incidence.shape[1],incidence.shape[0]
    t=nc*nr
    e=incidence.sum()
    if N==-1:
        if exact:
            N=ceil((e*(1-e/t)) *log((1-e/t)/accuracy) /2  )
        else:
            N=ceil((e/(2-2*e/t)) *log((1-e/t)/accuracy) )  
    if exact:
        return   N,<int>rewire_bipartite_ex(<unsigned short*>np.PyArray_DATA(incidence), nc,  nr, N, verbose, N*MAXITER)
    else:
         return  N,<int>rewire_bipartite(<unsigned short*>np.PyArray_DATA(incidence), nc,  nr, N, verbose)

def c_rewire_undirected(np.ndarray incidence,N=-1, verbose=1,  MAXITER=10, accuracy=0.00005,exact=True):
    cdef size_t e,nc,nr,t,d
    nc,nr= incidence.shape[1],incidence.shape[0]
    t=nc*nr/2
    e=incidence.sum()/2
    d=e/t
    if N==-1:
        if exact:
            N=ceil((e*(1-d)) *log((1-d)/accuracy) /2  )
        else:
            N=(e/(2*d^3-6*d^2+2*d+2))*log((1-d)/accuracy)
    if exact:
        return  N,<int>rewire_ex(<unsigned short*>np.PyArray_DATA(incidence), nc,  nr, N, verbose, N*MAXITER)
    else:
         return N,<int>rewire(<unsigned short*>np.PyArray_DATA(incidence), nc,  nr, N, verbose)

def c_analysis_bipartite(np.ndarray incidence,N=-1, verbose=1,  MAXITER=10, accuracy=0.00005,exact=True,step=10):
    cdef size_t e,nc,nr,t,dim
    cdef np.ndarray scores
    nc,nr= incidence.shape[1],incidence.shape[0]
    t=nc*nr
    e=incidence.sum()
    if N==-1:
        if exact:
            N=ceil((e*(1-e/t)) *log((1-e/t)/accuracy) /2  )
        else:
            N=ceil((e/(2-2*e/t)) *log((1-e/t)/accuracy) )  
    score=np.ascontiguousarray(np.zeros(N+1),dtype=np.double)
    if exact:
        dim=  <int>analysis_ex(<unsigned short*>np.PyArray_DATA(incidence), nc,nr,<double*>np.PyArray_DATA(score), step, N, verbose, MAXITER*N)
    else:
        dim=  <int>analysis(<unsigned short*>np.PyArray_DATA(incidence), nc,nr,<double*>np.PyArray_DATA(score), step, N, verbose)

    return N,score[0:dim]

def c_analysis_undirected(np.ndarray incidence,N=-1, verbose=1,  MAXITER=10, accuracy=0.00005,exact=True,step=10):
    cdef size_t e,nc,nr,t,dim,d
    cdef np.ndarray scores
    nc,nr= incidence.shape[1],incidence.shape[0]
    t=nc*nr/2
    e=incidence.sum()/2
    d=e/t
    if N==-1:
        if exact:
            N=ceil((e*(1-d)) *log((1-d)/accuracy) /2  )
        else:
            N=(e/(2*d^3-6*d^2+2*d+2))*log((1-d)/accuracy)
    score=np.ascontiguousarray(np.zeros(N+1),dtype=np.double)
    if exact:
        dim=  <int>analysis_undirected_ex(<unsigned short*>np.PyArray_DATA(incidence), nc,nr,<double*>np.PyArray_DATA(score), step, N, verbose, MAXITER*N)
    else:
        dim=  <int>analysis_undirected(<unsigned short*>np.PyArray_DATA(incidence), nc,nr,<double*>np.PyArray_DATA(score), step, N, verbose)

    return N,score[0:dim]

def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([("", a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

def incidence(np.ndarray x):
    edge_list=np.transpose((x==1).nonzero())
    edge_list[:,1]=edge_list[:,1]+x.shape[0]
    return i.Graph(list(edge_list))
class Rewiring:
    #data=None
    #N=None
    #data_rewired=None
    #jaccard_index=None
    #verbose=None
    #MAXITER=None
    #accuracy=None
    #step=None
    #exact=None
    #__type_of_data=None
    #__type_of_graph=None
    #__type_of_array=None
    def __init__(self,data,type_of_array=None,type_of_graph=None):
        self.data=data
        if type(self.data)==i.Graph:
            self.__type_of_data="graph"
            self.__type_of_array="None"
            if self.data.is_directed():
                print "Directed graph are not supported.\n" 
                self.data=None
            else:    
                if type_of_graph in ["bipartite","undirected"] :  
                    self.__type_of_graph=type_of_graph
                else:
                    print "The type of graph must be bipartite or undirected.\n" 
                    self.data=None
            if self.data.has_multiple():
                print "Multiple edge are not supported, I will simplify the graph. \n" 
                self.data=self.data.simplify()
        else:    
            if type(self.data)==np.ndarray:
                self.__type_of_data="array"
                if type_of_array not in ["edgelist_b","incidence","adjacence","edgelist_u"]:
                    print "The input type of array is not supported or must be given. \n" 
                    self.data=None
                else:    
                    self.__type_of_array=type_of_array
                    if self.__type_of_array in ["edgelist_b","incidence"]:
                        self.__type_of_graph="bipartite"
                    else:
                        self.__type_of_graph="undirected"
                if self.__type_of_array=="edgelist_b" or self.__type_of_array=="edgelist_u" :
                    if self.data.shape[1]!=2:
                        print "Edgelist must contain 2 colums" 
                        self.data=None
                        self.__type_of_array=None
                    else:
                        self.data=unique_rows(self.data)
                        self.data=np.ascontiguousarray(self.data,dtype=np.uintp)

                else:
                    self.data=np.ascontiguousarray(self.data,dtype="H")
            else:
                print "Data type not supported.\n" 
                self.data=None
        print "Object created: array="+self.__type_of_array+" data="+self.__type_of_data+" graph="+self.__type_of_graph     
    def rewire(self,N=-1,verbose=1,MAXITER=10, accuracy=0.00005,exact=True):
        self.N=N
        self.verbose=verbose
        self.MAXITER=MAXITER
        self.accuracy=accuracy
        self.exact=exact
        if self.__type_of_data=="graph":
            if self.__type_of_graph=="bipartite":
                ##get edgelist
                result=np.array(self.data.get_edgelist())
                left=np.ascontiguousarray(result[:,0],dtype=np.uintp)
                right=np.ascontiguousarray(result[:,1]-min(result[:,1]),dtype=np.uintp)
                tmp=c_rewire_sparse_bipartite(left,right,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact)
                result=np.vstack((left,right)).T
                result=i.Graph(list(result))
            else:
                #edgelist=np.ascontiguousarray(np.array(self.data.get_edgelist()),dtype=np.uintp)
                #tmp=c_rewire_undireced_sparse(edgelist,)                
                print "Not yet wrapped\n"
                return "Rewiring algorithm not performed"
        if self.__type_of_data=="array":
            if self.__type_of_array=="edgelist_u":
                #edgelist=np.copy(self.data)
                #tmp=c_rewire_undirected_sparse(edgelist,)
                print("Not yet wrapped\n")
                return "Rewiring algorithm not performed"
            if self.__type_of_array=="edgelist_b":
                result=np.copy(self.data)
                left=np.ascontiguousarray(result[:,0],dtype=np.uintp)
                right=np.ascontiguousarray(result[:,1]-min(result[:,1]),dtype=np.uintp)
                tmp=c_rewire_sparse_bipartite(left,right,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact)
                result=np.vstack((left,right)).T
                result=i.Graph(list(result))
            if self.__type_of_array=="incidence":
                result=np.ascontiguousarray(np.copy(self.data),dtype="H")
                tmp=c_rewire_bipartite(result,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact)
              
            if self.__type_of_array=="adjacence":
                result=np.ascontiguousarray(np.copy(self.data),dtype="H")
                tmp=c_rewire_undirected(result,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact)
        self.N=tmp[0]
        self.data_rewired=result
        if tmp[1]==0:
            return "Successfully rewired"
        else:
            return  "Rewiring algorithm not performed"
    def similarity(self):
            if self.data_rewired is None:
                print "First rewire the graph :-)."
                return -1
            else:
                if self.__type_of_array=="edgelist_b" or self.__type_of_array=="edgelist_u":
                    j_i=len(set(list(self.data)).intersection(set(list(self.data_rewired))))
                    return j_i/(len(list(self.data))-j_i)
                if self.__type_of_data=="graph":
                    m1=self.data.get_edgelist
                    m2=self.data_rewired.get_edgelist
                    j_i=len(set(m1).intersection(set(m2)))
                    return j_i/(len(list(m1))-j_i)
                if self.__type_of_array=="incidence" or self.__type_of_array=="adjacence":
                    return   (self.data*self.data_rewired).sum()/((self.data+self.data_rewired).sum()-(self.data*self.data_rewired).sum())

    def analysis(self,N=-1,verbose=1,MAXITER=10, accuracy=0.00005,exact=True,step=10): 
        #i can conmpute efficently the ji only incidence/adjacency matrix
        self.verbose=verbose
        self.MAXITER=MAXITER
        self.accuracy=accuracy
        self.exact=exact 
        self.step=step
        self.N=N
        result=np.ascontiguousarray(np.copy(self.data),dtype="H")

        if self.__type_of_graph=="bipartite" and self.__type_of_array=="incidence":
            tmp=c_analysis_bipartite(result,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact,self.step)        
        else:  
            if self.__type_of_graph=="unidrected" and self.__type_of_array=="adjacence":
                tmp=c_analysis(result,self.N, self.verbose,  self.MAXITER, self.accuracy,self.exact,self.step)      
            else:
                print "Give me an incidence or adjacency matrix"    
                return -1
        self.data_rewired=result
        self.jaccard_index=tmp[1]
        self.N=tmp[0]
        return "Analysis completed."
    def sampler(self,path,K=2000,max=1000,N=-1,verbose=0,MAXITER=10, accuracy=0.00005,exact=True):
        if not os.path.exists(path):
            os.makedirs(path)
        num_sub=ceil(K/max)
        for i in range(0,num_sub):
            if not os.path.exists(path+"/"+str(i)):
                os.makedirs(path+"/"+str(i))
            for j in range(0,max):
                
                self.rewire(N=N,verbose=verbose,MAXITER=MAXITER, accuracy=accuracy,exact=exact)
                self.data=self.data_rewired
                out_file = path+"/"+str(i)+"/"+str(i)+"_"+str(j)
                if self.__type_of_data=="array":
                    np.save(arr=self.data_rewired,file=out_file)
                else:
                    np.save(arr=np.array(self.data_rewired.get_edgelist),file=out_file)
                #out_file.write(self.data_rewired)
                #out_file.close()
            print "Saved "+str(max*(i+1))+" files\n"
        print "Finished."
def read_BRCA(file):
    cr = csv.reader(open(file))
    lista=[]
    i=0
    rownames=[]
    for row in cr:
        if i !=0:
            lista.append(row[1:])
            rownames.append(row[0])
        else:
            colnames=row[1:]
        i=i+1    
    return np.array(lista),colnames,rownames