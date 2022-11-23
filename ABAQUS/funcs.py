#-*- coding:utf-8 -*-
import numpy as np
from .struct import Object
from odbAccess import *
import time
import os
##读取点集与单元集的函数
def get_input(file):
    P=Object()
    f=open(file,'r')
    A=f.readlines()
    node_list=[]
    ele_list=[]
    section_list=[]
    material_list=[]
    for (i,a) in enumerate(A):
        if '*End Part' in a:
            k=i
        if '*Material' in a:
            material_list.append(i)
    for (i,a) in enumerate(A[:k]):
        if '*Nset' in a:
            node_list.append(i)
        if '*Elset' in a:
            ele_list.append(i)
        if '** Section' in a:
            section_list.append(i)
    #添加点集
    for n in node_list:
        nname=A[n][:-1].split(',')[1].split('=')[-1]
        if 'generate' in A[n]:
            a=A[n+1][:-1].split(',')
            amin=int(a[0])
            amax=int(a[1])
            nodes=np.arange(amin,amax+1,1)
            P.add_nset(nname, nodes, 1)
        else:
            nodes=[]
            i=n+1
            while 1:
                a=A[i][:-1]
                if '*' in a:
                    break
                for am in a.split(','):
                    if am=='':
                        continue
                    nodes.append(int(am))
                i=i+1
            nodes=np.array(nodes)
            P.add_nset(nname, nodes, 0)
    #添加单元集
    for e in ele_list:
        elname=A[e][:-1].split(',')[1].split('=')[-1]
        if 'generate' in A[e]:
            a=A[e+1][:-1].split(',')
            amin=int(a[0])
            amax=int(a[1])
            elements=np.arange(amin,amax+1,1)
            P.add_elset(elname, elements, 1)
        else:
            elements=[]
            i=e+1
            while 1:
                a=A[i][:-1]
                if '*' in a:
                    break
                for am in a.split(','):
                    if am=='':
                        continue
                    elements.append(int(am))
                i=i+1
            elements=np.array(elements)
            P.add_elset(elname, elements, 0)
    #添加区域
    for s in section_list:
        secname=A[s][:-1].split(': ')[-1]
        elname=A[s+1][:-1].split(',')[1].split('=')[-1]
        mname=A[s+1][:-1].split(',')[-1].split('=')[-1]
        P.add_section(secname,mname,elname)
    #添加材料
    for m in material_list:
        mname=A[m][:-1].split(',')[-1].split('=')[-1]
        elastic=float(A[m+4][:-1].split(',')[0])
        stress=[5e6,1e7,2e7]
        P.add_material(mname, elastic, stress,0)
    f.close()
    return P

def modify_input(P,file1,file2,ELASTIC,F):
    f1=open(file1,'r')
    f2=open(file2,'w')
    A=f1.readlines()
    for (i,a) in enumerate(A):
        if '*End Part' in a:
            k1=i
        if '** MATERIALS' in a:
            k2=i
        if '** PREDEFINED FIELDS' in a:
            k3=i
        if '** Name: shear1   Type: Surface traction' in a:
            k4=i
    for (i,a) in enumerate(A[:k1]):
        if '*Nset, nset=C1' in a:
            k0=i
    for i in range(k0):
        f2.write(A[i])
    for k in P.nset.keys():       
        if P.ngen[k]:
            a='*Nset, nset={:s},generate\n'.format(k)
            amin=np.min(P.nset[k])
            amax=np.max(P.nset[k])
            a=a+'{:d},{:d}, 1\n'.format(amin,amax)
            f2.write(a)
        else:
            a='*Nset, nset={:s}\n'.format(k)
            nodes=P.nset[k]
            lis=np.arange(0,len(nodes),16)
            lis=list(lis)
            lis.append(len(nodes))
            for l in range(len(lis)-1):
                nods=nodes[lis[l]:lis[l+1]]
                a=a+','.join('%5d'%n for n in nods)+'\n'
            f2.write(a)
    for k in P.elset.keys():        
        if P.elgen[k]:
            a='*Elset, elset={:s}, generate\n'.format(k)
            amin=np.min(P.elset[k])
            amax=np.max(P.elset[k])
            a=a+'{:d},{:d}, 1\n'.format(amin,amax)
            f2.write(a)
        else:
            a='*Elset, elset={:s}\n'.format(k)
            elements=P.elset[k]
            lis=np.arange(0,len(elements),16)
            lis=list(lis)
            lis.append(len(elements))
            for l in range(len(lis)-1):
                eles=elements[lis[l]:lis[l+1]]
                a=a+','.join('%5d'%n for n in eles)+'\n'
            f2.write(a)
    for k in P.sections.keys():
        a='** Section: {:s}\n*Solid Section, elset={:s}, material={:s}\n,\n'.\
            format(k,P.sections[k].elset_name,P.sections[k].material_name)
        f2.write(a)
    for i in range(k1,k2+2):
        f2.write(A[i])
    #改变材料信息
    Elastic1=ELASTIC[0]
    Elastic2=ELASTIC[1]
    for i in range(Elastic1.shape[0]):
        name='Inv1_{:d}'.format(i)
        P.change_material(name,Elastic1[i])
    for i in range(Elastic2.shape[0]):
        name='Inv2_{:d}'.format(i)
        P.change_material(name,Elastic2[i])
    #
    #写入材料信息
    for k in P.materials.keys():
        if not P.materials[k].pla:
            a='*Material, name={:s}\n*Density\n2500.,\n*Elastic\n'.format(k)
            a=a+'{:8e}, 0.25\n'.format(P.materials[k].elastic)
        else:
            stress=P.materials[k].stress
            a='*Material, name={:s}\n*Density\n2500.,\n*Elastic\n'.format(k)
            a=a+'{:8e}, 0.25\n'.format(P.materials[k].elastic)
            a=a+'*Mohr Coulomb\n20., 0.1\n*Mohr Coulomb Hardening\n'
            a=a+'{:8e},0.\n{:8e},0.0002\n{:8e},0.001\n'.format(stress[0],stress[1],stress[2])
        f2.write(a)
    #写入切应力信息
    for i in range(k3,k4):
        f2.write(A[i])
    a='** Name: shear1   Type: Surface traction\n*Dsload, op=NEW\n'
    a=a+'Y, TRSHR, {:8e}, 1., 0., 0.\n'.format(F[0])
    a=a+'** Name: shear2   Type: Surface traction\n*Dsload, op=NEW\n'
    a=a+'Y, TRSHR, {:8e},0., 0., 1.\n'.format(F[1])
    f2.write(a)
    for i in range(k4+6,len(A)):
        f2.write(A[i])
    f1.close()
    f2.close()
    
def read_out(file,Mutuli_Node):
    o=openOdb(file)
    class Node:
        def __init__(self,c1,c2):
            self.c1=c1   #前后两个点
            self.c2=c2
            self.distance=np.sqrt(sum((c1-c2)**2))
            self.U=[]
        def original_U(self,U0):
            self.U0=U0
        def add_U(self,u):
            self.U.append(u)
    Coord={}   #表示每个点的坐标值
    nodes=o.rootAssembly.instances['SOIL-1'].nodes
    node_label=[]
    for n in nodes:
        node_label.append(n.label)
    for k in Mutuli_Node.keys():
        c1=nodes[node_label.index(Mutuli_Node[k][0])].coordinates
        c2=nodes[node_label.index(Mutuli_Node[k][1])].coordinates
        Coord[k]=Node(c1,c2)   ##获取初始坐标值
    steps=o.steps.keys()
    for step in range(len(steps)):
        vall=o.steps[steps[step]].frames[-1].fieldOutputs['U'].values
        vall_label=[]
        for v in vall:
            vall_label.append(v.nodeLabel)
        for k in Mutuli_Node.keys():
            if step==Mutuli_Node[k][-2]:
                c1=Coord[k].c1
                c2=Coord[k].c2
                l1=vall[vall_label.index(Mutuli_Node[k][0])].data
                l2=vall[vall_label.index(Mutuli_Node[k][1])].data
                u=np.sqrt(sum((c1+l1-c2-l2)**2))
                Coord[k].original_U(u)
            if step in Mutuli_Node[k][-1]:
                l1=vall[vall_label.index(Mutuli_Node[k][0])].data
                l2=vall[vall_label.index(Mutuli_Node[k][1])].data
                c1=Coord[k].c1
                c2=Coord[k].c2
                u=np.sqrt(sum((c1+l1-c2-l2)**2))
                u=u-Coord[k].U0
                Coord[k].add_U(1000*u)
    Sim_U={}
    for k in Coord.keys():
        Sim_U[k]=np.array(Coord[k].U)
    del Coord
    o.close()
    return Sim_U

def finished(path):
    file=path+'/test.cid'
    while 1:
        if os.path.exists(file):
            time.sleep(1)
        else:
            break
    time.sleep(1)

def ab_num(x):
    if x>0:
        y=1
    elif x==0:
        y=0
    else:
        y=-1
    return y

def Out_puts(real_U,Sim_U,Iter,path,parallel):
    #parallel 表示是否并行
    import matlab.engine
    eng=matlab.engine.start_matlab()
    ERROR={}
    Error=0
    for k in real_U.keys():
        u0=real_U[k]
        u1=Sim_U[k]
        error=np.sum(abs((u1-u0)/u0))/u0.shape[0]
        error=error*100
        ERROR[k]=error
        Error=Error+error
    p=bool(parallel)
    if p:
        return Error 
     
    if Iter==0:
        f=open('Error.txt','w')
        a='Origin:\n'
    else:
        f=open('Error.txt','a')
        a='The {:d} Iter:\n'.format(Iter)
    for k in ERROR.keys():
        a=a+k+': {:4f}'.format(ERROR[k])+'\n'
    f.write(a)
    f.close()
    
    f0=open('plot/out.dat','w')
    if Iter==0:
        f=open('U.txt','w')
        a='Origin:\n'
    else:
        f=open('U.txt','a')
        a='The {:d} Iter:\n'.format(Iter)
    am=''
    for k in Sim_U.keys():
        a=a+k+':'+' '.join('%5f' %m for m in Sim_U[k])+'\n'
        am=am+k+':'+' '.join('%5f' %m for m in Sim_U[k])+'\n'
    f.write(a)
    f.close()
    f0.write(am)
    f0.close()
    eng.plot_error(path)
    return ERROR

def get_Q(crack_index,l):
    n=crack_index.shape[0]
    Q=np.zeros((n,n))
    for i in range(Q.shape[0]):
        for j in range(Q.shape[1]):
            a1=crack_index[i,-1]
            a2=crack_index[j,-1]
            d=a1-a2
            Q[i,j]=np.exp(-d/l)
    return Q
    