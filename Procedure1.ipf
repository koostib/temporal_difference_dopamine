//Stetsenko & Koos, Code 1 supplement
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <Function Grapher>
#include "Macintosh HD:Users:jaime:Library:Preferences:WaveMetrics:Igor Pro 6 Intel:Packages:WMFunctionGrapher:MyFunc"

// NOTES TO USER
// Initialize with running 'doit()' after loading the procedure .  Code is for forward Euler method with 0.001s time steps
// To run simulation of cell responses to V(t) square inputs of increasing amplitude execute 'doit()'.  Output is the 'out' array (see bellow) for cell activities etc..  
// The weighted sum continuous  output tis in the 'da_w_1_n' waves where n is the trial number. 
// V(t) inputs  are generated in  'doit()'  as the the wave 'input2';  the code is set up to generate square input pulses 
// To generate the stochastic synaptic  output to DA neurons execute 'discr(n)'  where the input variable n determines the number of trials. 
// The waves 't4_n' and 't5_n' are the stochastic synaptic outputwaves, where n is the serial number of the trial 
// Multiple trials are concatenated in waves 'wdest1' and 'wdest2' and when prompted by the save popup window, save  as't4.dat' and 't5.dat'.  The created files are text files with appropriate header and format for use as is in XPPAUT 
// Notes for function 'myint()'	
	//Most variables are  defined in the function 'doit()'
	// 'yw' array contains the variable values updated at each step
	// 'dydx' array are the momentary differential values for each variable
	// 'out' is a 2D array that contains each variable value for each time point
	// variable 0 is the normalized activity of T2Ns; Variable 1 is the T3N activity
	// Variable 2 is the synaptic efficacy of the facilitating synapse 
	// Variable 3-5 are the GABAB    simulation variables
	// Variable 6 is the synaptic effciacy of the depressing synapse
	// Initial values are commentized bellow; These  need to be decomentized when the first trial is run, after that the steady state end of each trial  can be used to initialize the run, if the same lines are commentized.


// Function to integrate the ODEs using the forward Euler method describing the cells' responses to different inputs. dt is variable step_size [s] and is set in function 'doit()' (0.1 or 1 ms was used for cell responses; 1ms for stochastic)
Function myInt(run)
variable run
NVAR n_sites
wave dydx
wave yw
variable i, sg
NVAR k_R_activation
NVAR k_R_deactivation
NVAR k_R_recovery
NVAR k_R_desense
NVAR k_G_activate
NVAR k_G_deactivate
NVAR Kd_GABAb
NVAR step_size
NVAR k_S_deplete
NVAR k_S_recovery
NVAR sim_sel
variable TR
variable GABA_CArry
wave input, out, GABAB
wave input2, input, input3
wave intensity, vals, mem
Wave pw	
string cwn0, cwn1
make /O/N=8 dydx

dydx=0
out=0
GABAB=0.439
dydx[3]=0
dydx[4]=0
dydx[5]=0
GABAB=(1/0.25)*yw[5]^n_sites/(Kd_GABAb+yw[5]^n_sites)
GABA_CArry=GABAB
out=0
sg=-1
yw[6]=1-0.885

// Initial conditions; uncomment for first run
// yw[0]=0.557089
// yw[1]=0.72455
// yw[2]=0.163351
// yw[3]=0.025657
// yw[4]=0.0359301
// yw[5]=0.016038
// yw[6]=0.5// 0.490815
// yw[7]=0.238189

For (i=0; i<(numpnts(GABAB)); i+=1)

	dydx[0] =(step_size)*(1000/pw[0])*(1.2*pw[2] -yw[0]- 1*(GABA_CArry))   
	dydx[1] =(step_size)*(1000/pw[1])*(0.9*pw[3]-yw[2]*0.18*yw[0]-yw[1]  -input2[i] )  
	dydx[2]=(step_size)*(1000/pw[4])*((facil(yw[0]))-yw[2]   )*0.8
	dydx[3]=(step_size)* ( (1- yw[3]- yw[4])*k_R_activation* (yw[1]) -yw[3]*k_R_deactivation  +k_R_recovery *yw[4]) 
	dydx[4]=(step_size)* ( k_R_desense* yw[3]  - k_R_recovery * yw[4] )  
	dydx[5]=(step_size)*(   k_G_activate* yw[3] -  k_G_deactivate* yw[5] )  
	dydx[6]=(step_size)*( (1- yw[6])*k_S_deplete* yw[1]-k_S_recovery* yw[6])

	yw[0]+=dydx[0]
	yw[1]+=dydx[1]
	yw[2]+=dydx[2]
	yw[3]+=dydx[3]
	yw[4]+=dydx[4]
	yw[5]+=dydx[5]
	yw[6]+=dydx[6]
	yw[7]+=dydx[7]

	out[i][0]=yw[0]
	out[i][1]=yw[1]
	out[i][2]=yw[2]
	out[i][3]=yw[3]
	out[i][4]=yw[4]
	out[i][5]=yw[5]

	if (sim_sel==6)
		out[i][6]=(1-yw[6])*yw[1]
	else
		out[i][6]=(1-yw[6])
	EndIf
	out[i][7]=yw[7]
	mem[i]=0.17*yw[0]
	if (Kd_GABAb!=-(yw[5]^n_sites))
		GABAB[i]=(1/0.25)*yw[5]^n_sites/(Kd_GABAb+yw[5]^n_sites)
	EndIf
	GABA_CArry=GABAB[i]
EndFOr
Doupdate
cwn0="cw_0_"+num2str(run)
cwn1="cw_1_"+num2str(run)
Make /O/N=(numpnts(GABAB)) $cwn0
Make /O/N=(numpnts(GABAB)) $cwn1
wave wcp = $cwn0
For (i=0; i<(numpnts(GABAB)); i+=1)
	wcp[i]=out[i][0]
EndFor
wave wcp = $cwn1
For (i=0; i<(numpnts(GABAB)); i+=1)
	wcp[i]=out[i][1]
EndFor
SetScale /P x, 0, 0.001, "s", $cwn0
SetScale /P x, 0, 0.001, "s", $cwn1

End

//Function to calculate steady state synaptic efficacy of the facilitating synapse as a function of presynaptic activity.   
function facil(rate)
variable rate
variable w
variable h=0.6
variable s=0.028

w=0.25*(1/(1+exp((-rate+h)/s))) 
if (rate<0.0)
	w=0.00
else
	w=(0.9/(1+exp((-rate+h)/s)))+0.1
endif
return w
End

// Function to generate the stochastic outputs for n trials, using 'spikeTrainExp2()' and 'spikeTrainExp()'
Function discr(n)
variable n
NVAR step_size
variable i
wave t4, t5
string wn=""
string wn1=""
string wn2=""
wave wdest1, wdest2
Make /O/N=(5/step_size) t4, t5
for (i=0; i<n; i+=1)
	 wn1="t4_"+num2str(i)
	Duplicate /O t4, $wn1
	 wn2="t5_"+num2str(i)
	 Duplicate /O t4, $wn2
	spikeTrainExp(wn1)
	spikeTrainExp2(wn2)
EndFor
lstr(n)
wdest1[0]=numpnts(wdest1)
wdest1[1]=0
wdest1[2]=numpnts(wdest1)
wdest2[0]=numpnts(wdest2)
wdest2[1]=0
wdest2[2]=numpnts(wdest2)
Save/G/M="\n" wdest1 as "t4.dat"
Save/G/M="\n" wdest2 as "t5.dat"
End


//  Next 2 functions  generate stochastic spike trains and IPSC conductance trains from the continuous activities of the T2 and T3Ns computed in 'myInt()'
// T3N output
Function spikeTrainExp(wns)
string wns
variable i, j, a, a1, a2, f, k1,u, k, s, dt, T,integ, p, interval_rand , shift, b, count, last, roll, cum_int, spk_count, th, gamma_scale, gamma_shape, isi, cont_syn, int_sum, noise,sp_count, n_released
NVAR release_p
NVAR step_size
NVAR  length
NVAR sim_sel
NVAR k_S_recovery
wave IPSPs, da_w_1_18, IPSP_single, out, IPSPs2, t4, spk_INT, spk, st
Make /O/N=(5/step_size) $wns
Make /O/N=(15/step_size) IPSPs, IPSPs2, IPSP_single

IPSPs=0
IPSPs2=0
IPSP_single=0
integ=0
Sim_sel=6
doit()
wave w = $wns
 
For (i=0; i<(length/step_size); i+=1)
	IPSP_single[i]=-(exp(-i/1.5) - exp(-i/15))* 2.01
EndFor

Make /O/N=1000 st
st=0


 
Print "T3N synaptic out"
For (k=0; k<45; k+=1)
print k, "  of 45"
integ=0
last=0
cum_int=0
spk_count=0
interval_rand=0
int_sum=0
st=0
sp_count=0

For (i=0; i<(length/step_size); i+=1)
	f=15*out[i][1]/0.66
	int_sum=(f)/1000
	if ((enoise(0.5)+0.5) < (int_sum))
		st[sp_count]=i
		sp_count+=1
		int_sum=0
	endif
EndFor

for (i=1; i<sp_count; i+=1)
	if ((st[i]-st[i-1])<4)
		st[i]=st[i]+4
	endif
endfor

sp_count=0
count=60
n_released=0
cont_syn=60*0.24
For (i=interval_rand; i<(length/step_size); i+=1)
	if (i>st[sp_count])
		sp_count+=1
		n_released=0
		count=round(cont_syn)
		for (b=0; b<round(cont_syn); b+=1)
			if ((enoise(0.5)+0.5) < (release_p))
				if (count>=1)
				count-=1
				n_released+=1
				EndIf
			EndIf
		EndFor //b 		
		cont_syn=	count
		For (j=i;j<(length/step_size); j+=1)
			if (((j)>=0) && ((j) < (i+250)) )
				IPSPs[j]+=IPSP_single[j-i]*(n_released)
			EndIf
		EndFor //j
	endif
	cont_syn+=step_size*k_S_recovery*(60-cont_syn)
EndFor
 EndFor
t4=0
For (i=0; i<5000; i+=1)
w[i]=IPSPs[i+7000]
EndFor
End


Function spikeTrainExp2(wns)
string wns
variable i, j, a, a1, a2, f, k1,u, k, s, dt, T,integ, p, interval_rand , shift, b, count, last, roll, cum_int, spk_count, th, gamma_scale, gamma_shape, isi, cont_syn, int_sum, noise,sp_count, n_released
NVAR release_p
NVAR step_size
NVAR  length
NVAR sim_sel
NVAR k_S_recovery
wave IPSPs, da_w_1_22, IPSP_single, out, IPSPs2, t4, spk_INT, spk, st
Make /O/N=5000 $wns
Make /O/N=(15/step_size) IPSPs, IPSPs2, IPSP_single
IPSPs=0
IPSPs2=0
IPSP_single=0
integ=0
Sim_sel=6
doit()
wave w = $wns
For (i=0; i<(length/step_size); i+=1)
	IPSP_single[i]=-(exp(-i/1.5) - exp(-i/15))*2.01
EndFor

Make /O/N=1000 st
st=0


Print "T2N synaptic out" 
For (k=0; k<35; k+=1)
print k, "  of 35"
integ=0
last=0
cum_int=0
spk_count=0
interval_rand=0
int_sum=0
st=0
sp_count=0
For (i=0; i<(length/step_size); i+=1)
	f=15*out[i][0]/0.54
	int_sum=(f)/1000
	if ((enoise(0.5)+0.5) < (int_sum))
		st[sp_count]=i
		sp_count+=1
		int_sum=0
	endif
EndFor
 

for (i=1; i<sp_count; i+=1)
if ((st[i]-st[i-1])<4)
st[i]=st[i]+4

endif
endfor


sp_count=0
count=5
n_released=0
cont_syn=5*0.24
For (i=interval_rand; i<(length/step_size); i+=1)
	if (i>st[sp_count])
		sp_count+=1
		n_released=0
		count=round(cont_syn)
		for (b=0; b<round(cont_syn); b+=1)
			if ((enoise(0.5)+0.5) < (0.5))
				if (count>=1)
					count-=1
					n_released+=1
				EndIf
			EndIf
		EndFor //b 		
		cont_syn=	count
		For (j=i;j<(length/step_size); j+=1)
			if (((j)>=i) && ((j) < (i+250)) )
				IPSPs[j]+=IPSP_single[j-i]*(n_released)
			EndIf
		EndFor //j
	endif
	cont_syn+=step_size*200*k_S_recovery*(5-cont_syn)
EndFor
EndFor
t4=0
For (i=0; i<5000; i+=1)
w[i]=IPSPs[i+7000]
EndFor
End



// Generates  the NMDA conductance waveform
Function makeNMDA(n)
variable n
NVAR step_size
variable i, j, dis_tau
dis_tau=140
wave NMDA
NMDA=0
for (i=0; i<n; i+=1)
for (j=0; j<1000; j+=1)
	NMDA[j+i*(5/step_size)+(3/step_size)]=(exp(-1*(j)/dis_tau)-exp(-1*(j)/5))
EndFor
EndFor
NMDA[0]=numpnts(NMDA)
NMDA[1]=0
NMDA[2]=numpnts(NMDA)
End

(5/step_size)

// Concatenates the individual trial outputs into a single long array
function lstr(n)
variable n
string /G s_name1=""
string /G s_name2=""
variable i
string wn1
string wn2

wave wdest1
for (i=0; i<n; i+=1)
 wn1="t4_"+num2str(i)+";"
s_name1+=wn1
endfor

concatenate /NP/O s_name1, wdest1
for (i=0; i<n; i+=1)
	wn2="t5_"+num2str(i)+";"
	s_name2+=wn2
endfor
print s_name2
concatenate /NP/O s_name2, wdest2
End


// Function for generating the response trajectories for phase portraits
Function traj()
NVAR step_size
NVAR length
variable i
wave trajw_1,  trajw_2, out
Make/O/N=(length/step_size) trajw_1,  trajw_2
for (i=(0); i<(length/step_size); i+=1)
	trajw_1[i]=out[i][0]//*0.22
	 trajw_2[i]=out[i][1]
EndFor
End


//**************************************  MAIN  ************************************** 


Function doit()
variable i, j
wave pw	
variable   tau1=20
variable  tau2=20
variable /G ratemax1=1.0
variable /G ratemax2=0.75
variable /G k_recovery=1
variable /G length=15
variable /G sim_sel=6 
variable /G n_sites=1.3
variable /G k_R_activation=12/4
variable /G k_R_deactivation=360/4
variable /G k_R_recovery=30/4
variable /G k_G_activate=30/4
variable /G k_G_deactivate=48/4
variable /G k_R_desense=42/4
variable /G Kd_GABAb=(5e-2)^n_sites
variable /G step_size=0.001 // [s]
variable /G  k_S_deplete=3.75
variable /G  k_S_recovery=0.65
variable /G release_p=1-e^(-k_S_deplete*0.6581/15)
variable slope1=0.5
variable slope2=0.5
variable slowFB=800
 variable slowSSn, cnt
wave input, input3
wave errW, wave2, wb
wave ips
wave c1_fix
wave c2_fix
wave intensity, linear
wave vals
wave GABAB
wave xw
wave yw
wave out
wave diff , mem
string wave_name=""
string wave_name1=""
string wave_name2=""
wave vals
wave GABAB
wave xw
wave yw
wave out, shit
wave diff , mem

make /O/N=12 intensity
make /O/N=2 errW
errW[0]=0.00001
errW[1]=0.00001

Make /O/N=24 ips
print 4/step_size
make /O/N=(length/step_size , 8) out,  copy_out
make /O/N=(length/step_size) input, diff, c1_fix, c2_fix, GABAB, mem, input3
make /O/N=12 pw
make /O/N=12 vals
make/O/N=13 wave2, wb
make/O/N=21 linear 
vals=0
ips=0
pw=0
make /O/N=(length/step_size) input2
SetScale /P x 0, step_size,"s", input, input2, c1_fix, c2_fix, GABAB, mem, input3
input=0
SetScale /P x 0, step_size, "s", out,copy_out, diff
SetScale /P x 0.3, 0.2, "s", wave2
input2=0
make /O/N=8  yw

pw[0]= tau1
pw[1]= tau2
pw[2] =ratemax1
pw[3]= ratemax2
pw[4]= slowFB
pw[5]= 155
pw[6]=slowFB
yw[0]=0.6
yw[1]=0.3945
yw[2]=0.174
yw[3]=0.7
yw[4]=0.183
cnt=0
For (j=1; j<20; j+=3)
input=0
input2=0

// V(t) inputs, optogenetic currents etc. are generated here.  The V(t) input is 'input2'; The code is for 2s long square pulse inputs between 2s to 10s 

for (i=(8/step_size); i<(10/step_size); i+=1)  // delay: for (i=(1/step_size); i<((1.1+(j-1)*0.1)/step_size); i+=1)
	input2[i]=0.0188406* j 
EndFor

input=0
for (i=(8.5/step_size); i<(10.0/step_size); i+=1)
input[i]=j*2/(45*3.6)
EndFor

input3=0
for (i=(9.0/step_size); i<(10.5/step_size); i+=1)
	input3[i]=   1   //  exponential input3[i]=   (1.46*2.6*6/(45*4.6))  *   exp(-1*(i-3/step_size) /  ( 0.3/step_size)   )  
EndFor

wave_name1 = "da_w_1_"+num2str(j)
wave_name2 = "da_w_2_"+num2str(j)
myInt(j)

for (i=0; i<(length/step_size); i+=1)
	diff[i]=0.12*out[i][0]+1*out[i][6]
EndFor
ips[j]=1-(diff[    ((1.146+(j-1)*0.1)/step_size)]-0.57183)/0.57183

Duplicate /O diff, $wave_name1

traj()
EndFor
End


