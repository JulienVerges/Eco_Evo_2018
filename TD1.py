import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

styles=["solid",(0,(3,5,1,5,1,5)),(0, (3, 10, 1, 10, 1, 10)), (0, (3, 5, 1, 5)), (0, (5, 5)), (0, (1, 5))]

R,C = sp.symbols(("R","C"))
r,K,a,Th,q,d=sp.symbols(("r","K","a","Th","q","d"))
#r, Th, q, K, d, a
standard = (0.2,1,1,100,0.25,0.02)

"""Tracé des isoclines et du champ d(R,C)/dt """
def isoclines(R_number,C_number, R_max=100, C_max=20, param=standard):
	r_,Th_,q_,K_,d_,a_ = param
	R_range=np.linspace(0,R_max,R_number)
	C_range=np.linspace(0,C_max,C_number)
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	f1=expr1[0].subs({r:r_, Th:Th_,q:q_,K:K_,d:d_,a:a_})
	f2=expr2[0].subs({r:r_, Th:Th_,q:q_,K:K_,d:d_,a:a_})
	Rnullcline_func = sp.lambdify(R,f1,dummify=0)
	Cnullcline_func = sp.lambdify(C,f2,dummify=0)
	plt.plot(R_range,Rnullcline_func(R_range))
	plt.plot([Cnullcline_func(0) for i in C_range],C_range)

def champ_vecteurs(R_number,C_number, R_max=100, C_max=20, param=standard):
	r_,Th_,q_,K_,d_,a_ = param	
	R_range=np.linspace(0,R_max,R_number)
	C_range=np.linspace(0,C_max,C_number)
	xx,yy= np.meshgrid(R_range, C_range)
	u,v= r_*(1-xx/K_) - a_*xx*yy/(1+a_*Th_*xx) , q_*a_*xx*yy/(1+a_*Th_*xx) - d_*yy
	plt.quiver(xx,yy,u,v)

"""Calcul approché d'une solution du système
dt : pas temporel
skip : nombres de points à oublier au début de la dynamique"""
def dynamique_pop(R_0,C_0,dt=0.1,param=standard,n_iter=10000,skip=9000):
	R=R_0
	C=C_0
	t=[k*dt for k in range(n_iter)]
	def my_system(var,t,r,Th,q,K,d,a):
		R,C=var
		return([r*R*(1-R/K)-a*R*C/(1+a*Th*R) , q*a*R*C/(1+a*Th*R)-d*C])
	sol=odeint(my_system, [R_0,C_0], t, param)
	solution=sol[skip::]
	x=[u[0] for u in solution]
	y=[u[1] for u in solution]
	return(x,y)

"""Effets de la variation des paramètres sur le comportement asymptotique
Les variables ind_style, ind_couleur, n_styles ont une utilité purement graphique"""
def variation_r(R_0,C_0,r_min,r_max,r_number):
	couleurs1=["navy","blue","dodgerblue"]
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	r_values=np.linspace(r_min,r_max,r_number)
	print(r_values)
	ligne=[]
	for r in r_values:
		param=(r,standard[1],standard[2],standard[3],standard[4],standard[5])
		x,y=dynamique_pop(R_0,C_0,0.1,param)
		l=plt.plot(x,y,label="r="+str(r),color=couleurs1[ind_couleur])
		ligne.append(l)
		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.legend()
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	

def variation_K(R_0,C_0,K_min,K_max,K_number):
	couleurs1=["navy","blue","dodgerblue"]
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	K_values=np.linspace(K_min,K_max,K_number)
	ligne=[]
	for K in K_values:
		param=(standard[0],standard[1],standard[2],K,standard[4],standard[5])
		x,y=dynamique_pop(R_0,C_0,0.1,param)
		l=plt.plot(x,y,label="K="+str(K),color=couleurs1[ind_couleur],linestyle=styles[ind_style])
		ligne.append(l)
		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.legend()
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	

def variation_d(R_0,C_0,d_min,d_max,d_number):
	couleurs1=["navy","blue","dodgerblue"]
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	d_values=np.linspace(d_min,d_max,d_number)
	print(d_values)
	ligne=[]
	for d in d_values:
		param=(standard[0],standard[1],standard[2],standard[3],d,standard[5])
		x,y=dynamique_pop(R_0,C_0,0.1,param)
		l=plt.plot(x,y,label="d="+str(d),color=couleurs1[ind_couleur],linestyle=styles[ind_style])
		ligne.append(l)
		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.legend()
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	

def variation_Th(R_0,C_0,Th_min,Th_max,Th_number):
	couleurs1=["navy","blue","dodgerblue"]
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	Th_values=np.linspace(Th_min,Th_max,Th_number)
	ligne=[]
	for Th in Th_values:
		param=(standard[0],Th,standard[2],standard[3],standard[4],standard[5])
		x,y=dynamique_pop(R_0,C_0,0.1,param)
		l=plt.plot(x,y,label="Th="+str(Th),color=couleurs1[ind_couleur],linestyle=styles[ind_style])
		ligne.append(l)
		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.legend()
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	

def variation_a(R_0,C_0,a_min,a_max,a_number):
	couleurs1=["navy","blue","dodgerblue"]
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	a_values=np.linspace(a_min,a_max,a_number)
	print(a_values)
	ligne=[]
	for a in a_values:
		param=(standard[0],standard[1],standard[2],standard[3],standard[4],a)
		x,y=dynamique_pop(R_0,C_0,0.1,param)
		l=plt.plot(x,y,label="a="+str(a),color=couleurs1[ind_couleur],linestyle=styles[ind_style])
		ligne.append(l)
		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.legend()
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	

"""Effets de la variation des paramètres sur les isoclines"""
def isoclines_r(r_min,r_max,r_number):
	couleurs1=["navy","blue","dodgerblue"]
	couleurs2=["firebrick","red","orangered"]
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	R_range=np.linspace(0,100,0.01)
	C_range=np.linspace(0,80,0.01)
	r_range=np.linspace(r_min,r_max, r_number)
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	for r_var in r_range:
		f1=expr1[0].subs({r:r_var, Th:1,q:1,K:100,d:0.25,a:0.02})
		f2=expr2[0].subs({r:r_var, Th:1,q:1,K:100,d:0.25,a:0.02})
		Rnullcline_func = sp.lambdify(R,f1,dummify=0)
		Cnullcline_func = sp.lambdify(C,f2,dummify=0)
		x=list(R_range)
		y=list(Rnullcline_func(R_range))
		plt.plot(x,y,color=couleurs1[ind_couleur],linestyle=styles[ind_style],label = "dR/dt=0, r="+str(r_var))
		x=[Cnullcline_func(0) for i in C_range] 
		y=list(C_range)
		plt.plot(x,y,color=couleurs2[ind_couleur],linestyle=styles[ind_style],label = "dC/dt=0, r="+str(r_var))

		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	plt.legend()
	


def isoclines_K(K_min,K_max,K_number):
	couleurs1=["navy","blue","dodgerblue"]
	couleurs2=["firebrick","red","orangered"]
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	R_range=np.linspace(0,100,1000)
	C_range=np.linspace(0,20,1000)
	K_range=np.linspace(K_min,K_max, K_number)
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	for K_var in K_range:
		f1=expr1[0].subs({r:0.2, Th:1,q:1,K:K_var,d:0.25,a:0.02})
		f2=expr2[0].subs({r:0.2, Th:1,q:1,K:K_var,d:0.25,a:0.02})
		Rnullcline_func = sp.lambdify(R,f1,dummify=0)
		Cnullcline_func = sp.lambdify(C,f2,dummify=0)
		x=list(R_range)
		y=list(Rnullcline_func(R_range))
		plt.plot(x,y,color=couleurs1[ind_couleur],linestyle=styles[ind_style],label = "dR/dt=0, K="+str(K_var))
		x=[Cnullcline_func(0) for i in C_range] 
		y=list(C_range)
		plt.plot(x,y,color=couleurs2[ind_couleur],linestyle=styles[ind_style],label = "dC/dt=0, K="+str(K_var))

		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	plt.legend()

def isoclines_d(d_min,d_max,d_number):
	couleurs1=["navy","blue","dodgerblue"]
	couleurs2=["firebrick","red","orangered"]
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	R_range=np.linspace(0,100,1000)
	C_range=np.linspace(0,20,1000)
	d_range=np.linspace(d_min,d_max, d_number)
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	for d_var in d_range:
		f1=expr1[0].subs({r:0.2, Th:1,q:1,K:100,d:d_var,a:0.02})
		f2=expr2[0].subs({r:0.2, Th:1,q:1,K:100,d:d_var,a:0.02})
		Rnullcline_func = sp.lambdify(R,f1,dummify=0)
		Cnullcline_func = sp.lambdify(C,f2,dummify=0)
		x=list(R_range)
		y=list(Rnullcline_func(R_range))
		plt.plot(x,y,color=couleurs1[ind_couleur],linestyle=styles[ind_style],label = "dR/dt=0, d="+str(d_var))
		x=[Cnullcline_func(0) for i in C_range] 
		y=list(C_range)
		plt.plot(x,y,color=couleurs2[ind_couleur],linestyle=styles[ind_style],label = "dC/dt=0, d="+str(d_var))

		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	plt.legend()

def isoclines_Th(Th_min,Th_max,Th_number):
	couleurs1=["navy","blue","dodgerblue"]
	couleurs2=["firebrick","red","orangered"]
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	R_range=np.linspace(0,100,1000)
	C_range=np.linspace(0,20,1000)
	Th_range=np.linspace(Th_min,Th_max, Th_number)
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	for Th_var in Th_range:
		f1=expr1[0].subs({r:0.2, Th:Th_var,q:1,K:100,d:0.25,a:0.02})
		f2=expr2[0].subs({r:0.2, Th:Th_var,q:1,K:100,d:0.25,a:0.02})
		Rnullcline_func = sp.lambdify(R,f1,dummify=0)
		Cnullcline_func = sp.lambdify(C,f2,dummify=0)
		x=list(R_range)
		y=list(Rnullcline_func(R_range))
		plt.plot(x,y,color=couleurs1[ind_couleur],linestyle=styles[ind_style],label = "dR/dt=0, Th="+str(Th_var))
		x=[Cnullcline_func(0) for i in C_range] 
		y=list(C_range)
		plt.plot(x,y,color=couleurs2[ind_couleur],linestyle=styles[ind_style],label = "dC/dt=0, Th="+str(Th_var))

		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	plt.legend()

def isoclines_a(a_min,a_max,a_number):
	couleurs1=["navy","blue","dodgerblue"]
	couleurs2=["firebrick","red","orangered"]
	expr1=sp.solve(r*R*(1-R/K) - a*C*R/(1+a*Th*R),C)
	expr2 = sp.solve((q*a*R*C)/(1+a*Th*R)-d*C,R)
	R_range=np.linspace(0,100,1000)
	C_range=np.linspace(0,20,1000)
	a_range=np.linspace(a_min,a_max, a_number)
	n_styles=len(styles)
	ind_style=0
	ind_couleur=0
	for a_var in a_range:
		f1=expr1[0].subs({r:0.2, Th:1,q:1,K:100,d:0.25,a:a_var})
		f2=expr2[0].subs({r:0.2, Th:1,q:1,K:100,d:0.25,a:a_var})
		Rnullcline_func = sp.lambdify(R,f1,dummify=0)
		Cnullcline_func = sp.lambdify(C,f2,dummify=0)
		x=list(R_range)
		y=list(Rnullcline_func(R_range))
		plt.plot(x,y,color=couleurs1[ind_couleur],linestyle=styles[ind_style],label = "dR/dt=0, a="+str(a_var))
		x=[Cnullcline_func(0) for i in C_range] 
		y=list(C_range)
		plt.plot(x,y,color=couleurs2[ind_couleur],linestyle=styles[ind_style],label = "dC/dt=0, a="+str(a_var))

		if ind_style==n_styles-1:
			ind_style=0
			if ind_couleur==2:
				ind_couleur=0
			else:
				ind_couleur+=1
		else:
			ind_style+=1
	plt.xlabel("Ressources")
	plt.ylabel("Consommateur")
	plt.legend()





isoclines_a(0.01,0.03,15)
plt.show()
