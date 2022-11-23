#-*- coding:utf-8 -*-
import numpy as np
from odbAccess import *
import os
import shutil
import time
from ABAQUS.funcs import *
from abaqus import *
from abaqusConstants import *
import copy
##读取约束，并进行裂缝定位
Constraint=np.array([])
f=open('Constraint.txt','r')
A=f.readlines()
for a in A:
    am=[int(s) for s in a[:-1].split()]
    am=np.array(am).reshape(1,-1)
    if not Constraint.any():
        Constraint=am
    else:
        Constraint=np.concatenate((Constraint,am),axis=0)
f.close()
a=np.arange(Constraint.shape[0],0,-1)-1
Constraint=Constraint[a]   #倒序，从上到下
cm=[]
for i in range(Constraint.shape[0]):
    c=Constraint[i]
    if c[c<200].shape[0]!=0:
        cm.append(i)
cm=np.array(cm)
cm_index=np.min(cm)
for i in range(1,cm.shape[0]):
    if cm[i]-1 not in cm:
        crack_coun=cm[i]     #两条裂缝的分界线
        break
w=0
rock_index=[]
crack1_index=np.array([])
crack2_index=np.array([])
for i in range(Constraint.shape[0]):
    for j in range(Constraint.shape[1]):
        if Constraint[i,j]>200:
            rock_index.append(w)
        else:
            a=np.array([i,j,w,0]).reshape(1,-1)
            if i<crack_coun-5:
                if crack1_index.shape==(0,):
                    crack1_index=a
                else:
                    crack1_index=np.concatenate((crack1_index,a),axis=0)
            else:
                if crack2_index.shape==(0,):
                    crack2_index=a
                else:
                    crack2_index=np.concatenate((crack2_index,a),axis=0)
        w=w+1
#寻找周围像素点, crack的最后一维表示周围像素点的个数
a=Constraint.shape[0]
b=Constraint.shape[1]
for i in range(crack1_index.shape[0]):
    a=crack1_index[i,0]
    b=crack1_index[i,1]
    pixel=Constraint[a-2:a+3,b-2:b+3]
    crack1_index[i,-1]=pixel[pixel>200].shape[0]
for i in range(crack2_index.shape[0]):
    a=crack2_index[i,0]
    b=crack2_index[i,1]
    pixel=Constraint[a-2:a+3,b-2:b+3]
    crack2_index[i,-1]=pixel[pixel>200].shape[0]

##读取单元集合，进行定位  
x_min=-97.71
z_min=200.0
lx=34.77
lz=15.042
nx=80
nz=np.max(cm)-np.min(cm)+1
dx=lx/nx
dz=lz/nz
z_min=z_min-cm_index*dz

P=get_input('test.inp') 
o=openOdb('test.odb')
nodes=o.rootAssembly.instances['SOIL-1'].nodes
elements=o.rootAssembly.instances['SOIL-1'].elementSets['INVERSE'].elements
MG_elements=list(P.elset['MG'])
lay_elements=list(P.elset['lay4'])
node_label=[]
for n in nodes:
    node_label.append(n.label)
node_arg=np.argsort(node_label)
Center=np.zeros((len(elements),4))
i=0
for e in elements:
    ele=e.connectivity
    X=[]
    Z=[]
    for n in ele:
        X.append(nodes[node_arg[n-1]].coordinates[0])
        Z.append(nodes[node_arg[n-1]].coordinates[2])
    X=np.mean(X)
    Z=np.mean(Z)
    idx=(X-x_min)//dx
    idz=(Z-z_min)//dz
    if idx<0 or idx>nx or idz<0:
        l=[idx,idz,-1,e.label]
    else:
        ids=idz*nx+idx
        l=[idx,idz,ids,e.label]
    Center[i]=l
    i=i+1
#找出基岩与软弱夹层
stress=[3e6,1e7,2e7]
Crack_elements=[]
for i in range(crack1_index.shape[0]):
    l=crack1_index[i,-2]
    name='Inv1_{:d}'.format(i)
    element_crack=Center[:,-1][Center[:,-2]==l]
    P.add_elset(name,element_crack,0)
    P.add_material(name,0,stress,1)
    P.add_section(name,name,name)
    Crack_elements.extend(list(element_crack))   
for i in range(crack2_index.shape[0]):
    l=crack2_index[i,-2]
    name='Inv2_{:d}'.format(i)
    element_crack=Center[:,-1][Center[:,-2]==l]
    P.add_elset(name,element_crack,0)
    P.add_material(name,0,stress,1)
    P.add_section(name,name,name)
    Crack_elements.extend(list(element_crack))  
for e in Crack_elements:
    if e in MG_elements:
        MG_elements.remove(int(e))
    else:
        lay_elements.remove(int(e))
lay_elements=np.array(lay_elements)
MG_elements=np.array(MG_elements)
P.add_elset('lay4',lay_elements,0)
P.add_elset('MG',MG_elements,0)
o.close()

##真实条件
Object={}
f=open('Object.txt','r')
A=f.readlines()
for a in A:
    name=a[:-1].split(':')[0]
    am=a[:-1].split(':')[1].split()
    U=[]
    for an in am:
        U.append(float(an))
    Object[name]=np.array(U)
f.close()
Mutuli_Node={}
Mutuli_Node['M4008']=[101,74978,0,[1,2,3,4,5,6,7,8,10,11]]
Mutuli_Node['M4067']=[4614,72685,3,[4,5,6,7,8,10,11]]
Mutuli_Node['M4068']=[70,13898,3,[4,5,6,7,8,11]]

Section_Elastic={}
Section_Elastic['lay2']=2.7e10 
Section_Elastic['lay3']=2.6e10
Section_Elastic['lay4']=2.5e10
Section_Elastic['lay1']=(Section_Elastic['lay2']+Section_Elastic['lay3'])/2
Section_Elastic['lay5']=(Section_Elastic['lay4']+Section_Elastic['lay3'])/2
Section_Elastic['MG']=Section_Elastic['lay4']+1e9
for k in Section_Elastic.keys():
    P.change_material(k,Section_Elastic[k])
##运行初始正模型
if os.path.exists('plot'):
    shutil.rmtree('plot')
os.mkdir('plot')
pwd=os.getcwd()
if os.path.exists('Krun'):
    shutil.rmtree('Krun')
os.mkdir('Krun')
Elastic1=1.0e10*np.ones((crack1_index.shape[0],))
Elastic2=1.6e10*np.ones((crack2_index.shape[0],))
file1='test.inp'
file2='Krun/test.inp'
ELASTIC=[Elastic1,Elastic2]
F=np.array([0,0])
modify_input(P,file1,file2,ELASTIC,F)
os.chdir('Krun')
mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
mdb.jobs['test'].submit()
os.chdir(pwd)
del mdb.jobs['test']
time.sleep(15)
finished('Krun')
file='Krun/test.odb'
Sim_U=read_out(file,Mutuli_Node)
Error=Out_puts(Object,Sim_U,0,'Krun',1)
ERROR=Out_puts(Object,Sim_U,0,'Krun',0)
##反演循环
Iter_num=10
e_num=3    #反演的弹性模量数目
f_num=2     #反演的切应力数目
move_elastic=3e9
move_traction=8e4
Inv_elastic=['lay2','lay3','lay4']   ##需要反演的分区
Q1=get_Q(crack1_index,10)
[U1,S1,V1]=np.linalg.svd(Q1)
Q2=get_Q(crack2_index,10)
[U2,S2,V2]=np.linalg.svd(Q2)
u1=abs(U1[:,0])
u2=abs(U2[:,0])
delta=2e8/np.max(abs(u1))
var1=u1*delta
delta=2e8/np.max(abs(u2))
var2=u2*delta

for Iter in range(Iter_num):
    if os.path.exists('running'):
        shutil.rmtree('running')
    os.mkdir('running')
    PATH=[]
    ########运行正模型，求梯度矩阵
    for p in range(len(Inv_elastic)):
        path='running/'+Inv_elastic[p]
        PATH.append(path)
        os.mkdir(path)
        section_Elastic=Section_Elastic.copy()
        Q=copy.copy(P)
        section_Elastic[Inv_elastic[p]]=section_Elastic[Inv_elastic[p]]+2e8
        for k in section_Elastic.keys():
            Q.change_material(k,section_Elastic[k])
        file2=path+'/test.inp'
        modify_input(Q,file1,file2,ELASTIC,F)
        os.chdir(path)
        mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
        mdb.jobs['test'].submit()
        os.chdir(pwd)
        del mdb.jobs['test']
        del Q
    for p in range(f_num):
        path='running/traction'+str(p+1)
        PATH.append(path)
        os.mkdir(path)
        shear=F.copy()
        shear[p]=shear[p]+3e3
        file2=path+'/test.inp'
        Q=copy.copy(P)
        modify_input(Q,file1,file2,ELASTIC,shear)
        os.chdir(path)
        mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
                     numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
        mdb.jobs['test'].submit()
        os.chdir(pwd)
        del mdb.jobs['test']
        del Q
    ####调整两条裂缝
    elastic1=Elastic1+var1
    Ele=[elastic1,Elastic2]
    path='running/crack1'
    os.mkdir(path)
    PATH.append(path)
    file2=path+'/test.inp'
    Q=copy.copy(P)
    modify_input(Q,file1,file2,Ele,F)
    os.chdir(path)
    mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
             numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
    mdb.jobs['test'].submit()
    os.chdir(pwd)
    del mdb.jobs['test']
    del Q
    
    elastic2=Elastic2+var2
    Ele=[Elastic1,elastic2]
    path='running/crack2'
    os.mkdir(path)
    PATH.append(path)
    file2=path+'/test.inp'
    Q=copy.copy(P)
    modify_input(Q,file1,file2,Ele,F)
    os.chdir(path)
    mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
             numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
    mdb.jobs['test'].submit()
    os.chdir(pwd)
    del mdb.jobs['test']
    del Q
    for path in PATH:
        finished(path)
    
    #########调整反演参数
    Grad=[]
    for p in range(e_num):
        file='running/'+Inv_elastic[p]+'/test.odb'
        u=read_out(file,Mutuli_Node)
        error=Out_puts(Object,u,0,'Krun',1)
        grad=error-Error
        print(grad)
        Grad.append(error-Error)
    k=np.max([move_elastic*np.exp(-float(Iter)/3),5e8])/np.max(abs(np.array(Grad)))
    for p in range(e_num):
        Section_Elastic[Inv_elastic[p]]=Section_Elastic[Inv_elastic[p]]-Grad[p]*k
    Section_Elastic['lay1']=(Section_Elastic['lay2']+Section_Elastic['lay3'])/2
    Section_Elastic['lay5']=(Section_Elastic['lay4']+Section_Elastic['lay3'])/2
    Section_Elastic['MG']=Section_Elastic['lay4']+1e9
    for k in Section_Elastic.keys():
        P.change_material(k,Section_Elastic[k])
    #
    Grad=[]
    for p in range(f_num):
        file='running/'+'traction'+str(p+1)+'/test.odb'
        u=read_out(file,Mutuli_Node)
        error=Out_puts(Object,u,0,'Krun',1)
        Grad.append(error-Error)
    k=np.max([move_traction*np.exp(-float(Iter)/3),2e4])/np.max(abs(np.array(Grad)))
    for p in range(f_num):
        F[p]=F[p]-Grad[p]*k
    #
    Grad=[] 
    file='running/crack1/test.odb' 
    u=read_out(file,Mutuli_Node) 
    u0=Object['M4008']
    ua=u['M4008']
    error=100*np.sum(abs((ua-u0)/u0))/u0.shape[0]
    Error=ERROR['M4008']
    Grad.append(error-Error)
    file='running/crack2/test.odb' 
    u=read_out(file,Mutuli_Node)
    ua=u['M4008']      
    error=100*np.sum(abs((ua-u0)/u0))/u0.shape[0]
    Grad.append(error-Error)
    if Grad[0]>0:
        k=2e9/np.max(abs(u1))
        Elastic1=Elastic1-ab_num(Grad[0])*k*u1
        Elastic1[Elastic1<2e9]=2e9
    if Grad[1]>0:
        k=4e8/np.max(abs(u2))
        Elastic2=Elastic2-ab_num(Grad[1])*k*u2
    ELASTIC=[Elastic1,Elastic2]
        
    ###重新计算
    if os.path.exists('Krun'):
        shutil.rmtree('Krun')
    os.mkdir('Krun')
    file2='Krun/test.inp'
    modify_input(P,file1,file2,ELASTIC,F)
    os.chdir('Krun')
    mdb.JobFromInputFile(name='test',inputFileName='test.inp',parallelizationMethodExplicit=DOMAIN,
                         numDomains=8,multiprocessingMode=DEFAULT, numCpus=8)
    mdb.jobs['test'].submit()
    os.chdir(pwd)
    del mdb.jobs['test']
    time.sleep(15)
    finished('Krun')
    file='Krun/test.odb'
    Sim_U=read_out(file,Mutuli_Node)
    Error=Out_puts(Object,Sim_U,Iter+1,'Krun',1)
    ERROR=Out_puts(Object,Sim_U,Iter+1,'Iter{:d}'.format(Iter+1),0)
    
    