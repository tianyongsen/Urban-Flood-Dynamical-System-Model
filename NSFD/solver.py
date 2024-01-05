
import numpy as np
import math
import random
from numpy import linalg as LA

alpha=0.5 #衰减常量
b=10        #元胞边长，元胞为正方形
area=100
h=8
q=0.3
endtime=400
rownum=8
colnum=8

# Z0=[[random.randint(0,5) for j in range (colnum)] for i in range(rownum)]
# h0=[[random.randint(0,1) for j in range (colnum)] for i in range(rownum)]

Z0_original=((1,0.9,0.8,0.8,0.8,0.8,0.8,0.8),(0.9,0.7,0.7,0.7,0.7,0.7,0.7,0.7),
    (0.8,0.7,0.6,0.6,0.6,0.6,0.6,0.6),(0.8,0.7,0.6,0.5,0.5,0.5,0.5,0.5),(0.8,0.7,
    0.6,0.5,0.4,0.4,0.4,0.4),(0.8,0.7,0.6,0.5,0.4,0.3,0.3,0.3),(0.8,0.7,0.6,0.5,
    0.4,0.3,0.2,0.2),(0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.2))
h0_original=((0.5,0.5,0.3,0.3,0.3,0.3,0.3,0.3),(0.5,0.5,0.3,0.3,0.3,0.3,0.3,
    0.3),(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3),(0.3,0.3,0.3,0.3,0.3,0.3,
    0.3,0.3),(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3),(0.3,0.3,0.3,0.3,0.3,
    0.3,0.3,0.3),(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3),(0.3,0.3,0.3 ,0.3,
    0.3,0.3,0.3,0.3))
Z0_modefied=((1,0.9,0.8,0.8,0.8,0.8,0.8,0.8),(0.9,0.7,0.7,0.7,0.7,0.7,0.7,0.7),
    (0.8,0.7,0.6,0.6,0.6,0.6,0.6,0.6),(0.8,0.7,0.6,0.5,0.5,0.5,0.5,0.5),(0.8,0.7,
    0.6,0.5,0.4,0.4,0.4,0.4),(0.8,0.7,0.6,0.5,0.4,0.3,0.3,0.3),(0.8,0.7,0.6,0.5,
    0.4,0.3,0.2,0.2),(0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.2))   #将第一个元胞高程提升至1m

Z0=Z0_original
#Z0=Z0_modefied
h0=h0_original

def solver_Ht(t_index,h0,Z,A):
	#function:  给定元胞序列和时间序列，以及初始水深和水位，返回元胞序列每个时刻的水位
	#Ht:矩阵  行为元胞序列，列为时间序列  , h0 and Z and A:narray  ; h0:初始水深，Z：高程，A：连接矩阵
    c=alpha*b/area   # 常量    
    Ht=h0+Z   #初始水深       

    Lambda,Vectors=LA.eigh(A)             #求解A的特征值和特征向量  narray类型，已经归一化
                                           #其中eigh(A)是求实对称矩阵或者海森矩阵
                                            #eig(A)是求一般矩阵
                                            #特征值向量为：array([0.1,0]),降序排列
	#print(Lambda)
    #V_inv=LA.inv(Vectors)                #求解特征向量的逆 
    V_inv=np.transpose(Vectors)            #因为已经归一化且A对称，所以特征向量矩阵为正交矩阵
                                            #正交矩阵的逆就是其转置
    U0=V_inv@Ht   						 #得到U的初值条件 @为numpy矩阵一般乘法
    coe=U0            				     #待定系数就是 U0

    #先求出(0,0)元胞


	#求解Ut，   dU/dt=-c*lambda*U  ,,U_i=e^{-c \lambda_i t}
    Ut=np.zeros(shape=(len(coe),len(t_index)))      #声明Ut，未初始化，用来存储所有t的计算结果，每一列为每一个时刻
    for i in range(len(t_index)):
    	Ut[:,i]=coe*np.exp(-c*t_index[i]*Lambda)    

    #将要输出的元胞序号变换到Ht,得到t时刻的水位或水深'
    Ht=Vectors@Ut
    for i in range(len(t_index)):
        Ht[:,i]=Ht[:,i]-Z

    return Ht

def solver_Ht_piecewise(t_index,h0,Z,A):
    #function:  给定元胞序列和时间序列，以及初始水深和水位，返回元胞序列每个时刻的水位
    #Ht:矩阵  行为元胞序列，列为时间序列  , h0 and Z and A:narray  ; h0:初始水深，Z：高程，A：连接矩阵
    c=alpha*b/area   # 常量    
    Ht=h0+Z   #初始水深   
    
    #零点之前
    Lambda,Vectors=LA.eigh(A)     #求解A的特征值和特征向量
    print(max(Lambda))
    V_inv=np.transpose(Vectors)   #求解特征向量的逆 
    U0=V_inv@Ht                #得到U的初值条件
    coe=U0                               #待定系数就是 U0
    t_e=find_zero_point(Vectors[0,:],U0,Lambda,c,Z[0]) #找到使得(0,0)水深为0的时刻
    print(t_e)


    #求解Ut，   dU/dt=-c*lambda*U  ,,U_i=e^{-c \lambda_i t}
    Ht=np.zeros(shape=(len(coe),len(t_index)))      #声明Ht，未初始化，用来存储所有t的计算结果，每一列为每一个时刻
    for i in range(len(t_index)):
        if t_index[i]<t_e: #小于零点的时刻按照之前的计算
            Ht[:,i]=Vectors@(coe*np.exp(-c*t_index[i]*Lambda))-Z            #覆盖掉上面的声明 

    #===============================
    #零点之后
    H_mid=Vectors@(coe*np.exp(-c*t_e*Lambda)) #间断时刻的水位
    B=A[1:,1:]    #将原(0,0)元胞的行列截掉
    B[0][0]=2     #将与原(0,0)元胞相邻的元胞的邻域数量减1
    B[7][7]=2
    
    LambdaB,VB=LA.eigh(B)     #求解B的特征值和特征向量
    #print(LambdaB)
    VB_inv=np.transpose(VB)   #求解特征向量的逆 
    U_mid0=VB_inv@H_mid[1:]                #得到U的初值条件     
    for i in range(len(t_index)):
        if t_index[i]>t_e: #小于零点的时刻按照之前的计算
            Ht[0][i]=0                          
            #求得新解
            temp=VB@(U_mid0*np.exp(-c*(t_index[i]-t_e)*LambdaB)) 
            for j in range(len(temp)):
                Ht[j+1][i]=temp[j]-Z[j+1]           
                                                  
    return Ht
def find_zero_point(r1,u0,eignvalues,coe,z0):
    #r1:特征向量矩阵Q的第一行，u0=Q^T H_0 ,eignvalues:特征值向量（降序），coe：综合系数
    def f(t):   #定义函数f,找到函数f的零点
        return r1@(u0*np.exp(-coe*t*eignvalues))-z0 #@：向量内积，*：向量哈达玛积或者向量数乘
    #二分法找零点 取初始时间区间[0,400]，对于我的问题足够了
    a=0;b=400;eps=.00001
    print(f(a),f(b))
    while(b-a>eps):  #这一块很不严谨，但是基于我对这个问题的先验，足够解决我的问题了
        if(f(a)*f(b)<0):
            c=(a+b)/2.
            if f(c)>0:
                a=c
            elif f(c)<0:
                b=c
            else:
                return c
        else:
            print("时间区间是错的")

    return c     #找到零点

def NSFD_solver(t_index,cellfield,Ht):
    #暂时用最简单直接的方法计算这个特殊的案例
    #构建元胞域矩阵
    #对元胞域矩阵每个元胞遍历，计算其下一步水深
    Ht.clear()
    phi_h=(1-math.exp(-1*q*h))/q   #时间步长   (1-e^(-qh))/q

    for step in t_index:
        NSFD_reserve_result(cellfield,Ht)    #保存当前步长的结果 
        NSFD_advance_step(cellfield,phi_h)   #推进到下一步长


    # h_index=[n*h for n in range(int(t_index[-1]/h)+1)]    #一般设置h可以被各个t整除
    # k=0
    # for step in h_index:   # 循环次数
    #     if step==t_index[k]:
    #         NSFD_reserve_result(cellfield,Ht)
    #         k=k+1
    #     NSFD_advance_step(cellfield,phi_h)


def NSFD_advance_step(cf,phi_h):
    rownum=len(cf)
    colnum=len(cf[0])
    link_map={}               #用来存储g_ij 
    #遍历构建g_ij
    for i in range(rownum):
        for j in range(colnum):
            flux=0
            #判断元胞位置，元胞位置不同邻域不同
            if i==0 :
                if j==0:      
                    NSFD_cal_outer_gij(phi_h,link_map,cf,2,(i,j),(i,j+1),(i+1,j))
                elif j==colnum-1:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,2,(i,j),(i,j-1),(i+1,j))                   
                else:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,3,(i,j),(i+1,j),(i,j+1),(i,j-1))                   
            elif i==rownum-1:
                if j==0:        
                    NSFD_cal_outer_gij(phi_h,link_map,cf,2,(i,j),(i-1,j),(i,j+1))
                elif j==colnum-1:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,2,(i,j),(i-1,j),(i,j-1))                    
                else:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,3,(i,j),(i-1,j),(i,j+1),(i,j-1))                   

            else:  #平凡情况，全部有四邻域
                if j==0:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,3,(i,j),(i-1,j),(i+1,j),(i,j+1))
                elif j==colnum-1:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,3,(i,j),(i-1,j),(i+1,j),(i,j-1))                    
                else:
                    NSFD_cal_outer_gij(phi_h,link_map,cf,4,(i,j),(i-1,j),(i+1,j),(i,j-1),(i,j+1))
    #遍历更新
    for i in range(rownum):
        for j in range(colnum):
            if i==0:
                if j==0:
                    if (i,j,i,j+1) in link_map:   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j+1)]         #找到gij 
                    else:
                        gi1=-link_map[(i,j+1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i+1,j) in link_map :   
                        gi2=link_map[(i,j,i+1,j)]
                    else:
                        gi2=-link_map[(i+1,j,i,j)]
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2)    #更新该元胞状态
                elif j==colnum-1:
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]        #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i+1,j) in link_map :   
                        gi2=link_map[(i,j,i+1,j)]
                    else:
                        gi2=-link_map[(i+1,j,i,j)]
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2)    #更新该元胞状态
                else:
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]       #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i+1,j) in link_map :   
                        gi2=link_map[(i,j,i+1,j)]
                    else:
                        gi2=-link_map[(i+1,j,i,j)]
                    if (i,j,i,j+1) in link_map :   
                        gi3=link_map[(i,j,i,j+1)]
                    else:
                        gi3=-link_map[(i,j+1,i,j)]                    
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2+gi3)    #更新该元胞状态
            elif i==rownum-1:
                if j==0:
                    if (i,j,i,j+1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j+1)]         #找到gij 
                    else:
                        gi1=-link_map[(i,j+1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i-1,j) in link_map :   
                        gi2=link_map[(i,j,i-1,j)]
                    else:
                        gi2=-link_map[(i-1,j,i,j)]
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2)    #更新该元胞状态
                elif j==colnum-1:
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]         #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i-1,j) in link_map :   
                        gi2=link_map[(i,j,i-1,j)]
                    else:
                        gi2=-link_map[(i-1,j,i,j)]
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2)    #更新该元胞状态
                else:
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]        #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i-1,j) in link_map :   
                        gi2=link_map[(i,j,i-1,j)]
                    else:
                        gi2=-link_map[(i-1,j,i,j)]
                    if (i,j,i,j+1) in link_map :   
                        gi3=link_map[(i,j,i,j+1)]
                    else:
                        gi3=-link_map[(i,j+1,i,j)]                    
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2+gi3)    #更新该元胞状态
            else:
                if j==0:
                    if (i,j,i-1,j) in link_map :   
                        gi1=link_map[(i,j,i-1,j)]
                    else:
                        gi1=-link_map[(i-1,j,i,j)]
                    if (i,j,i,j+1) in link_map :   
                        gi2=link_map[(i,j,i,j+1)]
                    else:
                        gi2=-link_map[(i,j+1,i,j)]       
                    if (i,j,i+1,j) in link_map :   
                        gi3=link_map[(i,j,i+1,j)]
                    else:
                        gi3=-link_map[(i+1,j,i,j)]      
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2+gi3)    #更新该元胞状态
                elif j==colnum-1:
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]         #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)] 
                    if (i,j,i-1,j) in link_map :   
                        gi2=link_map[(i,j,i-1,j)]
                    else:
                        gi2=-link_map[(i-1,j,i,j)]
                    if (i,j,i+1,j) in link_map :   
                        gi3=link_map[(i,j,i+1,j)]
                    else:
                        gi3=-link_map[(i+1,j,i,j)]      
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2+gi3)    #更新该元胞状态
                else:                    
                    if (i,j,i,j-1) in link_map :   #如果键在字典dict里返回true，否则返回false
                        gi1=link_map[(i,j,i,j-1)]         #找到gij 
                    else:
                        gi1=-link_map[(i,j-1,i,j)]         #找到gij的转置，并取负号 变为正值
                    if (i,j,i-1,j) in link_map :   
                        gi2=link_map[(i,j,i-1,j)]
                    else:
                        gi2=-link_map[(i-1,j,i,j)]
                    if (i,j,i,j+1) in link_map :   
                        gi3=link_map[(i,j,i,j+1)]
                    else:
                        gi3=-link_map[(i,j+1,i,j)]       
                    if  (i,j,i+1,j) in link_map :   
                        gi4=link_map[(i,j,i+1,j)]
                    else:
                        gi4=-link_map[(i+1,j,i,j)]     
                    cf[i][j]=cf[i][j]+phi_h*(gi1+gi2+gi3+gi4)    #更新该元胞状态
def NSFD_cal_outer_gij(phi_h,link_map,cf,num,id0,id1,id2=(),id3=(),id4=()):

    if num==1: 
        ai1_=cal_out_flux(cf[id0[0]][id0[1]],cf[id1[0]][id1[1]],Z0[id0[0]][id0[1]],Z0[id1[0]][id1[1]])   
        if ai1<=0:
            link_map[(id0[0],id0[1],id1[0],id1[1])]=cf[id0[0]][id0[1]]*ai1_/(cf[id0[0]][id0[1]]-phi_h*(ai1_))
    elif num==2:
        ai1_=cal_out_flux(cf[id0[0]][id0[1]],cf[id1[0]][id1[1]],Z0[id0[0]][id0[1]],Z0[id1[0]][id1[1]])    
        ai2_=cal_out_flux(cf[id0[0]][id0[1]],cf[id2[0]][id2[1]],Z0[id0[0]][id0[1]],Z0[id2[0]][id2[1]])      
        outflux=0
        if ai1_<=0:   
            outflux=outflux+ai1_
        if ai2_<=0:
            outflux=outflux+ai2_
        if ai1_<=0:
            link_map[(id0[0],id0[1],id1[0],id1[1])]=cf[id0[0]][id0[1]]*ai1_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai2_<=0:
            link_map[(id0[0],id0[1],id2[0],id2[1])]=cf[id0[0]][id0[1]]*ai2_/(cf[id0[0]][id0[1]]-phi_h*outflux)
    elif num==3:
        ai1_=cal_out_flux(cf[id0[0]][id0[1]],cf[id1[0]][id1[1]],Z0[id0[0]][id0[1]],Z0[id1[0]][id1[1]])    
        ai2_=cal_out_flux(cf[id0[0]][id0[1]],cf[id2[0]][id2[1]],Z0[id0[0]][id0[1]],Z0[id2[0]][id2[1]])   
        ai3_=cal_out_flux(cf[id0[0]][id0[1]],cf[id3[0]][id3[1]],Z0[id0[0]][id0[1]],Z0[id3[0]][id3[1]])  
        outflux=0
        if ai1_<=0:   
            outflux=outflux+ai1_
        if ai2_<=0:
            outflux=outflux+ai2_
        if ai3_<=0:
            outflux=outflux+ai3_
        if ai1_<=0:
            link_map[(id0[0],id0[1],id1[0],id1[1])]=cf[id0[0]][id0[1]]*ai1_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai2_<=0:
            link_map[(id0[0],id0[1],id2[0],id2[1])]=cf[id0[0]][id0[1]]*ai2_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai3_<=0:
            link_map[(id0[0],id0[1],id3[0],id3[1])]=cf[id0[0]][id0[1]]*ai3_/(cf[id0[0]][id0[1]]-phi_h*outflux)

    else:
        ai1_=cal_out_flux(cf[id0[0]][id0[1]],cf[id1[0]][id1[1]],Z0[id0[0]][id0[1]],Z0[id1[0]][id1[1]])    
        ai2_=cal_out_flux(cf[id0[0]][id0[1]],cf[id2[0]][id2[1]],Z0[id0[0]][id0[1]],Z0[id2[0]][id2[1]])   
        ai3_=cal_out_flux(cf[id0[0]][id0[1]],cf[id3[0]][id3[1]],Z0[id0[0]][id0[1]],Z0[id3[0]][id3[1]])  
        ai4_=cal_out_flux(cf[id0[0]][id0[1]],cf[id4[0]][id4[1]],Z0[id0[0]][id0[1]],Z0[id4[0]][id4[1]])  
        outflux=0
        if ai1_<=0:   
            outflux=outflux+ai1_
        if ai2_<=0:
            outflux=outflux+ai2_
        if ai3_<=0:
            outflux=outflux+ai3_
        if ai4_<=0:
            outfulx=outflux+ai4_
        if ai1_<=0:
            link_map[(id0[0],id0[1],id1[0],id1[1])]=cf[id0[0]][id0[1]]*ai1_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai2_<=0:
            link_map[(id0[0],id0[1],id2[0],id2[1])]=cf[id0[0]][id0[1]]*ai2_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai3_<=0:
            link_map[(id0[0],id0[1],id3[0],id3[1])]=cf[id0[0]][id0[1]]*ai3_/(cf[id0[0]][id0[1]]-phi_h*outflux)
        if ai4_<=0:
            link_map[(id0[0],id0[1],id4[0],id4[1])]=cf[id0[0]][id0[1]]*ai4_/(cf[id0[0]][id0[1]]-phi_h*outflux)

def NSFD_reserve_result(cf,Ht):
    Htemp=[]
    for i in range (len(cf)):
        for j in range(len(cf[0])):
            Htemp.append(cf[i][j]/area)   #转为向量
    Ht.append(Htemp)                 #存储到整个大Ht矩阵中

def cal_out_flux(v_center,v2,Z_center,Z2):
    #从h_center的角度计算出流通量，若是入流则设为0单位m3/t    
    h_center=v_center/area
    h2=v2/area  
    H_center=h_center+Z_center
    H2=h2+Z2
    if(H_center<H2):
        return 1         #返回一个正值用于辨别
    flux=-1*b*alpha*(H_center-H2) #保持公式逻辑 不做简化 
    return flux
def construct_A(A,rownum,colnum):
    #function ：构建8*8元胞范围的连接矩阵
    #A:naarys               #
    for i in range(rownum):             #一行一行的循环 将（i,j）映射(i)*colnum+j,
        for j in range(colnum):          #（i,j）的邻域也按照这个规则映射 
            if i==0:   
                if j==0:       
                    A[(i)*colnum+j,(i)*colnum+j]=2  #对角线元素
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1     #link  下边
                    A[(i)*colnum+j,i*colnum+j+1]=-1   #link 右边
                elif j==colnum-1:
                    A[(i)*colnum+j,(i)*colnum+j]=2  #对角线元素
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1     #link 下方
                    A[(i)*colnum+j,i*colnum+j-1]=-1   #link  左边
                else:
                    A[(i)*colnum+j,(i)*colnum+j]=3  #对角线元素
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1     #link  下边
                    A[(i)*colnum+j,(i)*colnum+j+1]=-1   #link  右边                  
                    A[(i)*colnum+j,(i)*colnum+j-1]=-1   #link  左边               
            elif i==rownum-1:
                if j==0:        
                    A[(i)*colnum+j,(i)*colnum+j]=2  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link  上边
                    A[(i)*colnum+j,(i)*colnum+j+1]=-1   #link 右边
                elif j==colnum-1:
                    A[(i)*colnum+j,(i)*colnum+j]=2  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link 上方
                    A[(i)*colnum+j,(i)*colnum+j-1]=-1   #link  左边
                else:
                    A[(i)*colnum+j,(i)*colnum+j]=3  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link  上边
                    A[(i)*colnum+j,(i)*colnum+j+1]=-1   #link  右边                  
                    A[(i)*colnum+j,(i)*colnum+j-1]=-1   #link  左边
            else: 
                if j==0:
                    A[(i)*colnum+j,(i)*colnum+j]=3  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link  上边
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1   #link  下边                  
                    A[(i)*colnum+j,(i)*colnum+j+1]=-1   #link  右边
                elif j==colnum-1:
                    A[(i)*colnum+j,(i)*colnum+j]=3  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link  上边
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1   #link  下边                  
                    A[(i)*colnum+j,(i)*colnum+j-1]=-1   #link  左边
                else:                               #最平凡情况，在里面
                    A[(i)*colnum+j,(i)*colnum+j]=4  #对角线元素
                    A[(i)*colnum+j,(i-1)*colnum+j]=-1     #link  上边
                    A[(i)*colnum+j,(i+1)*colnum+j]=-1   #link  下边                  
                    A[(i)*colnum+j,(i)*colnum+j-1]=-1   #link  左边
                    A[(i)*colnum+j,(i)*colnum+j+1]=-1   #link  右边     
def construct_Cellfield(rownum,colnum,h0,Z):
    #cellfield:list    
    #function ：给定行列和初始值，创建相应的矩阵,并赋予初始值
    cellfield=[]
    for i in range(rownum):
        cellfield.append([0 for j in range(colnum)])   #创建cellfield
    #H=h0+Z     #得到水位  narray
    H=h0        #使用水深  
    for index in range(rownum*colnum): 
        i=int(index/colnum)
        j=index%colnum        
        cellfield[i][j]=H[index]*area


    return cellfield

def RungeKutta_solver(t_index,cellfield,Z_field,Ht):
    #这里使用4级4阶RungeKutta方法，但是在hi=0时做了一些校正,避免水深为负值
    #t_index: 报告时间节点 ，要求每个输入时间节点可以被h整除
    Ht.clear()
    dt=h       #时间步长，将全局变量拿下来
    caltime_index=[n*h for n in range(int(t_index[-1]/h)+1)] #计算时间节点
    k=0
    for time in caltime_index:   # 循环次数
        if time==t_index[k]:     #计算时与报告时相等，则报告结果
            RungeKutta_reserve_result(cellfield,Ht)
            k=k+1
        RungeKutta_advance_step(cellfield,Z_field,dt)  #每步推荐
def RungeKutta_advance_step(cellfield,Z_field,dt):
    #cellfield :记录水深 Z_field：记录高程
    c=b*alpha/area        #用全局变量给的 
    rownum=len(cellfield)
    colnum=len(cellfield[0])

    Hij=[[cellfield[i][j]+Z_field[i][j] for j in range(colnum)]for i in range(rownum)]


    #k1
    k1=[[0 for j in range(colnum)]  for i in range(rownum)]  #矩阵形式  
    RungeKutta_cal_k(k1,Hij,cellfield,dt,c)               #计算k1  
    #计算中间水位 
    RungeKutta_cal_intermediateHij(k1,Hij,cellfield,Z_field,dt/2.)   

    #k2 
    k2=[[0 for j in range(colnum)]  for i in range(rownum)]  #矩阵形式
    RungeKutta_cal_k(k2,Hij,cellfield,dt,c)               #计算k2  
    #计算中间水位 
    RungeKutta_cal_intermediateHij(k2,Hij,cellfield,Z_field,dt/2.)  

    #k3
    k3=[[0 for j in range(colnum)]  for i in range(rownum)]  #矩阵形式 
    RungeKutta_cal_k(k3,Hij,cellfield,dt,c)               #计算k3  
    #计算中间水位 
    RungeKutta_cal_intermediateHij(k3,Hij,cellfield,Z_field,dt)  

    #k4 
    k4=[[0 for j in range(colnum)]  for i in range(rownum)]  #矩阵形式 
    RungeKutta_cal_k(k4,Hij,cellfield,dt,c)               #计算k4

    #更新 
    K=[[k1[i][j]+2*k2[i][j]+2*k3[i][j]+k4[i][j] for j in range(colnum)]for i in range(rownum)]  #计算总K
    RungeKutta_cal_intermediateHij(K,Hij,cellfield,Z_field,dt/6)  #计算得到最终水位
    for i in range(len(Hij)):
        for j in range(len(Hij[0])): 
            cellfield[i][j]=Hij[i][j]-Z_field[i][j]   #更新水深值 

def RungeKutta_cal_flux_partial(Hc,hc,H1,h1):
    if Hc>H1 and hc<=0:
        return 0
    elif Hc<H1 and h1<=0:
        return 0
    else: 
        return Hc-H1

def RungeKutta_cal_k(k,Hij,cellfield,dt,c):
    rownum=len(cellfield)
    colnum=len(cellfield[0])
    for i in range(rownum):
        for j in range(colnum):            
            if i==0:
                if j==0:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])              
                elif j==colnum-1:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])              
                else:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])                
            elif i==rownum-1:
                if j==0:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])              
                elif j==colnum-1:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])              
                else:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])                
            else :
                if j==0:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])  
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])

                elif j==colnum-1:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])    
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])            
                else:
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j-1],\
                        cellfield[i][j-1])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i-1][j],\
                        cellfield[i-1][j])
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i][j+1],\
                        cellfield[i][j+1])  
                    k[i][j]=k[i][j]+RungeKutta_cal_flux_partial(Hij[i][j],cellfield[i][j],Hij[i+1][j],\
                        cellfield[i+1][j])
            #乘以系数
            k[i][j]=-c*k[i][j]







def RungeKutta_cal_intermediateHij(k,Hij,cellfield,Z_field,dt):
    for i in range(len(Hij)):
        for j in range(len(Hij[0])):
            Hij[i][j]=cellfield[i][j]+Z_field[i][j]+dt*k[i][j]       
            if Hij[i][j]<Z_field[i][j]:   #即预估水深小于0,则将预估水深置为0 
                Hij[i][j]=Z_field[i][j]                  #暂时不调节邻域接受到的通量
            
                


def RungeKutta_reserve_result(cellfield,Ht):
    Htemp=[cellfield[i][j] for i in range(len(cellfield)) for j in range(len(cellfield[0]))]            
    Ht.append(Htemp)                 #存储到整个大Ht矩阵中
