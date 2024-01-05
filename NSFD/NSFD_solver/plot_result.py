import matplotlib.pyplot as plt
import numpy as np

#结果绘图

def plot_result_three(Y1,Y2,Y3,x):
	fig, axs = plt.subplots()
	axs.set_ylim(-0.2, 0.8)
	axs.set_xlim(0,x[-1]+2)
	for i in range(len(Y1[:,0])): 
		if i==0:     #
			axs.plot(x, Y1[i,:],'r-',label='exact',linewidth=0.8)			
			axs.plot(x, Y3[i,:],'--b',label='RungeKutta',linewidth=1.0)
			axs.plot(x, Y2[i,:],'-.g',label='NSFD',linewidth=1.0)
			axs.annotate('cell(0,0)',fontsize=12,xy=(320,0.03),xytext=(320, 0.03))
			maxerror=0.
			print("the maxerror of NSFD:")
			for k in range(len(x)):
				a=abs(Y2[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a			
			print(maxerror)               #输出最大误差
			maxerror=0.
			print("the maxerror of RK:")
			for k in range(len(x)):
				a=abs(Y3[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a
			print(maxerror)				
		elif i==1:
			axs.plot(x, Y1[i,:],'r-',linewidth=0.8)			
			axs.plot(x, Y3[i,:],'--b',linewidth=1)
			axs.plot(x, Y2[i,:],'-.g',linewidth=1)
			axs.annotate('cell(1,1)',fontsize=12,xy=(320,0.23),xytext=(320, 0.23))
			maxerror=0.
			print("the maxerror of NSFD:")
			for k in range(len(x)):
				a=abs(Y2[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a			
			print(maxerror)               #输出最大误差
			maxerror=0.
			print("the maxerror of RK:")
			for k in range(len(x)):
				a=abs(Y3[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a
			print(maxerror)	
		elif i==2:
			axs.plot(x, Y1[i,:],'r-',linewidth=0.8)			
			axs.plot(x, Y3[i,:],'--b',linewidth=1)
			axs.plot(x, Y2[i,:],'-.g',linewidth=1)
			axs.annotate('cell(4,4)',fontsize=12,xy=(320,0.53),xytext=(320, 0.53))
			maxerror=0.
			print("the maxerror of NSFD:")
			for k in range(len(x)):
				a=abs(Y2[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a			
			print(maxerror)               #输出最大误差
			maxerror=0.
			print("the maxerror of RK:")
			for k in range(len(x)):
				a=abs(Y3[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a
			print(maxerror)	
		elif i==3: 
			axs.plot(x, Y1[i,:],'r-',linewidth=0.8)			
			axs.plot(x, Y3[i,:],'--b',linewidth=1)
			axs.plot(x, Y2[i,:],'-.g',linewidth=1)
			axs.annotate('cell(7,7)',fontsize=12,xy=(320,0.73),xytext=(320, 0.73))
			maxerror=0.
			print("the maxerror of NSFD:")
			for k in range(len(x)):
				a=abs(Y2[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a			
			print(maxerror)               #输出最大误差
			maxerror=0.
			print("the maxerror of RK:")
			for k in range(len(x)):
				a=abs(Y3[i,k]-Y1[i,k])
				if a>maxerror:
					maxerror=a
			print(maxerror)	
	axs.grid()
	axs.legend(fontsize=12)
	axs.set_xlabel('time (s)',fontsize=15)
	axs.set_ylabel('water depth (m)',fontsize=15)
	#axs.set_title('in')
	plt.savefig("../result_graph/first.jpg",dpi=2000)
	plt.show()
def plot_result_six(Y1,Y2,Y3,x,x2): 
	fig, axs = plt.subplots()
	axs.set_ylim(-0.2, 0.8)
	axs.set_xlim(0,x[-1]+2)
	#输入格式：axs.annotate('cell(0,0)',fontsize=12,xy=(320,0.03),xytext=(320, 0.03))
	c_a=[['cell(0,0)',12,(320,0.03),(320, 0.03)],
	# ['cell(1,1)',12,(320,0.23),(320, 0.23)],
	# ['cell(4,4)',12,(320,0.53),(320, 0.53)],
	['cell(1,1)',12,(320,0.26),(320, 0.23)],
	['cell(4,4)',12,(320,0.56),(320, 0.53)],
	['cell(7,7)',12,(320,0.73),(320, 0.73)]]	

	#先把第一个画出来，顺便画图例
	axs.plot(x, Y1[0,:],'r-',label='exact',linewidth=0.8)			
	axs.plot(x, Y3[0,:],'--b',label='RungeKutta',linewidth=1.0)
	axs.plot(x2, Y2[0,:],'-.g',label='NSFD',linewidth=1.0)
	axs.annotate(c_a[0][0],fontsize=c_a[0][1],xy=c_a[0][2],xytext=c_a[0][3])

	for i in range(1,len(Y1[:,0])):
		axs.plot(x, Y1[i,:],'r-',linewidth=0.8)			
		axs.plot(x, Y3[i,:],'--b',linewidth=1.0)
		axs.plot(x2, Y2[i,:],'-.g',linewidth=1.0)
		axs.annotate(c_a[i][0],fontsize=c_a[i][1],xy=c_a[i][2],xytext=c_a[i][3])

	axs.grid()
	axs.legend(fontsize=12)
	axs.set_xlabel('time (s)',fontsize=15)
	axs.set_ylabel('water depth (m)',fontsize=15)
	#axs.set_title('in')
	plt.savefig("first.jpg",dpi=2000)
	plt.show()



#从txt中读取文件并绘图
def read_file(M,file):
	beginning_flag="x"  #x_location  u rho p    

	with open(file,'r') as f:
		#定位到数据开始位置
		for line in f:
			if line[0]==beginning_flag:
				break
		#读取数据，保存到矩阵
		for line in f:			
			row=tuple([float(i) for i in line.split()])		
			M.append(row)
		f.close()








