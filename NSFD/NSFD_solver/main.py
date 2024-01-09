

import numpy as np
import solver
import plot_result
#ctrl+H 替换g     
def main():
	#初始值
	rownum=solver.rownum
	colnum=solver.colnum
	Z0_list=[solver.Z0[i][j] for i in range(rownum) for j in range (colnum)]
	
	h0_list=[solver.h0[i][j] for i in range(rownum) for j in range (colnum)]
	#nindex=[n for n in range(0,64)]
	#nindex=[0,9,18,27,36,45,54,63]
	nindex=[0,9,36,63]
	timeindex=[t for t in range(0,solver.endtime,solver.h)]
	phi_h=(1-np.exp(-solver.q*solver.h))/solver.q   #求解phi_h
	timeindex_NSFD=[t for t in np.arange(0,solver.endtime,phi_h)] #NSFD的时间点
	

	Z0=np.array(Z0_list)
	h0=np.array(h0_list)
	H=Z0+h0
	average=(np.sum(H)-Z0[0])/63.  #计算最终稳定水位 		   
	print(average)

	#求取精确解 
	A=np.zeros((rownum*colnum,rownum*colnum))      #创建连接矩阵
	solver.construct_A(A,rownum,colnum)         #填充连接矩阵
	#Ht=solver.solver_Ht(timeindex,h0,Z0,A)
	Ht=solver.solver_Ht_piecewise(timeindex,h0,Z0,A)

	A2=np.zeros((rownum*colnum,rownum*colnum))      #创建连接矩阵
	solver.construct_A(A2,rownum,colnum)         #填充连接矩阵
	Ht2=solver.solver_Ht_piecewise(timeindex_NSFD,h0,Z0,A2)


	#NSFD
	Ht_NSFD=[]       #每一行是一个时刻
	cellfield=solver.construct_Cellfield(rownum,colnum,h0,Z0)  #创建元胞域矩阵 
	solver.NSFD_solver(timeindex_NSFD,cellfield,Ht_NSFD)


	#RungeKutta 
	Ht_RungeKutta=[]	
	cellfield2=[ [solver.h0[i][j] for j in range(len(solver.h0[0]))] for i in range(len(solver.h0))]   #初始水深
	solver.RungeKutta_solver(timeindex,cellfield2,solver.Z0,Ht_RungeKutta)


	#绘图
	Ht_slice=np.empty(shape=(len(nindex),len(timeindex)))
	Ht2_slice=np.empty(shape=(len(nindex),len(timeindex_NSFD)))
	Ht_slice_NSFD=np.empty(shape=(len(nindex),len(timeindex_NSFD)))
	Ht_slice_RK=np.empty(shape=(len(nindex),len(timeindex)))
	for i in range(len(nindex)):
		Ht_slice[i,:]=Ht[nindex[i],:]    #对Ht做切片 
		Ht2_slice[i,:]=Ht2[nindex[i],:]  #对Ht2做切片
		for j in range(len(timeindex)):			
			Ht_slice_RK[i,j]=Ht_RungeKutta[j][nindex[i]]  #切片   Ht_NSFD为list,Ht_slice_NSFD为array
		for j in range(len(timeindex_NSFD)):
			Ht_slice_NSFD[i,j]=Ht_NSFD[j][nindex[i]]  #切片   Ht_NSFD为list,Ht_slice_NSFD为array

	plot_result.plot_result_six(Ht_slice,Ht_slice_NSFD,Ht_slice_RK,timeindex,timeindex_NSFD)
	

	#计算误差大小
	print("NSFD误差： RK误差") 
	for i in range(len(nindex)):
		e_NSFD=0.
		e_RK=0.
		for j in range(len(timeindex_NSFD)):
			e=abs(Ht2_slice[i,j]-Ht_slice_NSFD[i,j])
			if e_NSFD<e:
				e_NSFD=e

		for j in range(len(timeindex)):
			e=abs(Ht_slice[i,j]-Ht_slice_RK[i,j])
			if e_RK<e:
				e_RK=e
		print(e_NSFD,e_RK)



			











if __name__ == '__main__':
	main()



