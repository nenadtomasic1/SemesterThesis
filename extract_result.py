import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import sys
import array as arr
import scipy.stats as st


COST_l_HF=[28.8,230.4,1843.2,14745.6 ]
COST_l_LF=[1,8,64,512]
W_L_arr=np.asarray( COST_l_HF )/np.asarray( COST_l_LF)


def load_file(full_path_file_name, var_QoI):
	with open(full_path_file_name) as f:
		data = json.load(f)

	res=[]
	samples_dict=data['samples']
	
	for each in samples_dict:
			sample_obj=samples_dict[each]			
			obj_value=sample_obj['obj']

			if(  str( obj_value[var_QoI] ) !="inf" ):				
				res.append( float(obj_value[var_QoI]) )
    
	return res


#
# assume len(ANY_list)>>30
#
def split_into_train_and_test(ANY_list, slize):
	
	train_list_1=[]
	train_list_2=[]
	train_list_3=[]
	train_list_4=[]
	train_list_5=[]
	train_list_6=[]
	train_list_7=[]
	train_list_8=[]
	train_list_9=[]
	train_list_10=[]

	for i in range(len(ANY_list)):		
		train_list_1.insert(i, ANY_list[i][:slize])
		train_list_2.insert(i, ANY_list[i][slize:2*slize])
		train_list_3.insert(i, ANY_list[i][2*slize:3*slize])
		train_list_4.insert(i, ANY_list[i][3*slize:4*slize])
		train_list_5.insert(i, ANY_list[i][4*slize:5*slize])
		train_list_6.insert(i, ANY_list[i][5*slize:6*slize])
		train_list_7.insert(i, ANY_list[i][6*slize:7*slize])
		train_list_8.insert(i, ANY_list[i][7*slize:8*slize])
		train_list_9.insert(i, ANY_list[i][8*slize:9*slize])
		train_list_10.insert(i, ANY_list[i][9*slize:10*slize])

	return train_list_1,train_list_2,train_list_3,train_list_4,train_list_5,train_list_6,train_list_7,train_list_8,train_list_9,train_list_10



def get_r_l_star(rho_l, w_l):	
	return (-1 + math.sqrt( (w_l*rho_l*rho_l ) / (1 - (rho_l * rho_l ) ) ) )
	


def get_LAMBDA_l(rho_l, r_l):	
	return (1- ( r_l*math.pow(rho_l,2)/(1+r_l) ) )	


def get_Y_l_array(HF_arr):
	var_Y_l_HF_arr=[]
	for i in range(len(HF_arr)):
		if(i==0):
			var_Y_l_HF_arr.insert(i, np.mean( HF_arr[i]) )
		else:			
			var_Y_l_HF_arr.insert(i, np.mean(HF_arr[i]) - np.mean(HF_arr[i-1]))
	return var_Y_l_HF_arr




def get_variances_for_diffs_Y_l(values_arr):	
	result=[]
	for i in range(len(values_arr)):		
		if(i==0):
			result.insert( i, np.var( values_arr[i]) )
		else:
			len1=len(values_arr[i])
			len2=len(values_arr[i-1])	
			## we need this, due to the fact that some arrays have slightly less elements than others, in our case in the LF arrays have not the same number of elements
			## its slightly less precise, but acceptable 		
			result.insert( i, np.var( values_arr[i] ) + np.var( values_arr[i-1] )  - 2*np.cov( values_arr[i][:np.minimum(len1, len2)] , values_arr[i-1][:np.minimum(len1, len2)])[0][1] )			
	return result




def get_N_l_MLMF(EPSILON_VAL, l, VAR_Y_array, rho_l_array, cost_array_HF ):
	val=0
	for k in range(len(VAR_Y_array)):		
		r_k_star=get_r_l_star( rho_l_array[k], W_L_arr[k] )
		val= val + math.sqrt( (VAR_Y_array[k]*cost_array_HF[k])/(1 - math.pow( rho_l_array[l] ,2) ) )*get_LAMBDA_l( rho_l_array[k], r_k_star)		
	return (EPSILON_VAL)*val*( math.sqrt( (1- math.pow( rho_l_array[l] ,2) )*(VAR_Y_array[l]/cost_array_HF[l]) ) )





def get_EPSILON_val_MLMF(N_l, l, VAR_Y_array, rho_l_array, cost_array_HF):
	val=0
	for k in range(len(VAR_Y_array)):		
		r_k_star=get_r_l_star( rho_l_array[k], W_L_arr[k] )		
		val=val + math.sqrt( (VAR_Y_array[k]*cost_array_HF[k])/(1 - math.pow( rho_l_array[l] ,2) ) )*get_LAMBDA_l( rho_l_array[k], r_k_star)	
	return N_l/( val*( math.sqrt( (1- math.pow( rho_l_array[l] ,2) )*( VAR_Y_array[l]/cost_array_HF[l] ) ) ) )




def calculate_VAR_Q_MLMF(epsilon_MLMF_val, VAR_Y_array, rho_l_array, cost_array_HF, cost_array_LF):	
	val=0
	comp_cost_value=0
	N_VALUE_arr=[]	
	N_VALUE_arr_pure=[]
	r_VALUE_arr=[]
	lambda_arr=[]
	for k in range(len(VAR_Y_array)):		
		if k==0:
			N_l_mlmf_val=get_N_l_MLMF(epsilon_MLMF_val, k, VAR_Y_array, rho_l_array, cost_array_HF)
			N_VALUE_arr_pure.insert(k,N_l_mlmf_val)
		else:	
			N_l_mlmf_val=math.ceil( get_N_l_MLMF(epsilon_MLMF_val, k, VAR_Y_array, rho_l_array, cost_array_HF) )
			N_VALUE_arr_pure.insert(k,get_N_l_MLMF(epsilon_MLMF_val, k, VAR_Y_array, rho_l_array, cost_array_HF))
		
		N_VALUE_arr.insert(k,N_l_mlmf_val )
		r_k_star=get_r_l_star( rho_l_array[k], W_L_arr[k] )
		r_VALUE_arr.insert(k, r_k_star)			
		lambda_arr.insert(k,  get_LAMBDA_l( rho_l_array[k], r_k_star) )
		val= val + ( VAR_Y_array[k]/N_l_mlmf_val )* get_LAMBDA_l( rho_l_array[k], r_k_star)		
		comp_cost_value=comp_cost_value+N_l_mlmf_val*(cost_array_HF[k]+cost_array_LF[k]*(1+r_k_star))
	return val,comp_cost_value,N_VALUE_arr,r_VALUE_arr,N_VALUE_arr_pure,lambda_arr



def get_N_l_MLMC1(EPSILON_VAL, l, VAR_Y_array, cost_array):
	varianz_y_values=np.array(VAR_Y_array,dtype=np.float)
	cost_array_values=np.array(cost_array,dtype=np.float)	
	return (EPSILON_VAL)*sum( np.sqrt(varianz_y_values*cost_array_values) )*( math.sqrt( (VAR_Y_array[l]/(cost_array[l] ) ) ) )
	


def getMLMC_final_VAR(EPSILON_VAL,VAL_arr, cost_array):
	N_VALUE_arr=[]
	comp_cost_value=0
	variances_val_arr=get_variances_for_diffs_Y_l(VAL_arr)
	
	val=0.0

	for i in range(len( variances_val_arr )): 		
		if(i==0):			
			N_VALUE=int( round( get_N_l_MLMC1(EPSILON_VAL, i,variances_val_arr, cost_array) ) )	
		else:			
			N_VALUE=get_N_l_MLMC1(EPSILON_VAL, i,variances_val_arr, cost_array)		
		N_VALUE_arr.insert(i,N_VALUE )		
		comp_cost_value+=N_VALUE*cost_array[i]
		val=val+( variances_val_arr[i]/(  N_VALUE ) ) 
	
	return val,comp_cost_value,N_VALUE_arr






def get_RHO_l_s(HF_arr,LF_arr):
	t1_HF,t2_HF,t3_HF,t4_HF,t5_HF,t6_HF,t7_HF,t8_HF,t9_HF,t10_HF=split_into_train_and_test(HF,10)
	t1_LF,t2_LF,t3_LF,t4_LF,t5_LF,t6_LF,t7_LF,t8_LF,t9_LF,t10_LF=split_into_train_and_test(LF,49)

	t1_HF_y=get_Y_l_array(t1_HF)
	t2_HF_y=get_Y_l_array(t2_HF)
	t3_HF_y=get_Y_l_array(t3_HF)
	t4_HF_y=get_Y_l_array(t4_HF)
	t5_HF_y=get_Y_l_array(t5_HF)
	t6_HF_y=get_Y_l_array(t6_HF)
	t7_HF_y=get_Y_l_array(t7_HF)
	t8_HF_y=get_Y_l_array(t8_HF)
	t9_HF_y=get_Y_l_array(t9_HF)
	t10_HF_y=get_Y_l_array(t10_HF)

	t1_LF_y=get_Y_l_array(t1_LF)
	t2_LF_y=get_Y_l_array(t2_LF)
	t3_LF_y=get_Y_l_array(t3_LF)
	t4_LF_y=get_Y_l_array(t4_LF)
	t5_LF_y=get_Y_l_array(t5_LF)
	t6_LF_y=get_Y_l_array(t6_LF)
	t7_LF_y=get_Y_l_array(t7_LF)
	t8_LF_y=get_Y_l_array(t8_LF)
	t9_LF_y=get_Y_l_array(t9_LF)
	t10_LF_y=get_Y_l_array(t10_LF)

	res_val=[]
	for i in range(4):		
		corr_val=np.corrcoef( [ t1_HF_y[i],t2_HF_y[i],t3_HF_y[i],t4_HF_y[i],t5_HF_y[i],t6_HF_y[i],t7_HF_y[i],t8_HF_y[i],t9_HF_y[i],t10_HF_y[i]  ] ,[ t1_LF_y[i],t2_LF_y[i],t3_LF_y[i],t4_LF_y[i],t5_LF_y[i],t6_LF_y[i],t7_LF_y[i],t8_LF_y[i],t9_LF_y[i],t10_LF_y[i] ] )[0,1]
		res_val.insert(1, corr_val )
	return res_val





def load_all_files(a_val_QoI):	
	HF_l3_list=load_file('./HF-l3-IsoDAR60MeV_samples_0.json', a_val_QoI )
	LF_l3_list=load_file('./LF-200/smallDt-LF-l3-IsoDAR60MeV_samples_0.json', a_val_QoI )
	HF_l2_list=load_file('./HF-l2-IsoDAR60MeV_samples_0.json', a_val_QoI )
	LF_l2_list=load_file('./LF-200/smallDt-LF-l2-IsoDAR60MeV_samples_0.json', a_val_QoI )
	HF_l1_list=load_file('./HF-l1-IsoDAR60MeV_samples_0.json', a_val_QoI )
	LF_l1_list=load_file('./LF-200/smallDt-LF-l1-IsoDAR60MeV_samples_0.json', a_val_QoI )
	HF_l0_list=load_file('./HF-l0-IsoDAR60MeV_samples_0.json', a_val_QoI )
	LF_l0_list=load_file('./LF-200/smallDt-LF-l0-IsoDAR60MeV_samples_0.json', a_val_QoI )
	HF_arr=[HF_l0_list,HF_l1_list,HF_l2_list,HF_l3_list]
	LF_arr=[LF_l0_list,LF_l1_list,LF_l2_list,LF_l3_list]	
	return HF_arr,LF_arr 

############# MAIN #############

val_QoI='HALO_X' 


HF,LF=load_all_files(val_QoI)
CORRELATIONS_RHO_l=get_RHO_l_s(HF,LF)
EPSILON_MLMF=get_EPSILON_val_MLMF(  len(HF[0]) , 0, get_variances_for_diffs_Y_l(HF) , CORRELATIONS_RHO_l, COST_l_HF)



if( (len(sys.argv) == 3) and  (sys.argv[1]=='-e') ):
	EPSILON_MLMF=2/( float( sys.argv[2] )*float( sys.argv[2] ) )



### do the calculations
result_mlmc_var,comp_cost_MLMC,n_samples_mlmc_arr=getMLMC_final_VAR(EPSILON_MLMF,HF,COST_l_HF)
result_mlmf_var,comp_cost_MLMF,n_samples_mlmf_var,r_values_list,n_samples_mlmf_var_pure,lambda_aval_arr=calculate_VAR_Q_MLMF( EPSILON_MLMF, get_variances_for_diffs_Y_l(HF), CORRELATIONS_RHO_l, COST_l_HF, COST_l_LF )

## display the results

print("HF array " )
print(" mean  HF[0]: " + str( np.mean( HF[0]) ).rjust(20) + "  var HF[0]: " + str( np.var(HF[0])).rjust(20)  + "   len HF : " + str( len(HF[0])) ) 
print(" mean  HF[1]: " + str( np.mean( HF[1]) ).rjust(20) + "  var HF[1]: " + str( np.var(HF[1])).rjust(20)  + "   len HF : " + str( len(HF[1])) )
print(" mean  HF[2]: " + str( np.mean( HF[2]) ).rjust(20) + "  var HF[2]: " + str( np.var(HF[2])).rjust(20)  + "   len HF : " + str( len(HF[2])) )
print(" mean  HF[3]: " + str( np.mean( HF[3]) ).rjust(20) + "  var HF[3]: " + str( np.var(HF[3])).rjust(20)  + "   len HF : " + str( len(HF[3])) )
print("")
print("LF array " )
print(" mean LF[0]: " + str( np.mean( LF[0]) ).rjust(20) + "  var LF[0]: " +  str( np.var(LF[0])).rjust(20)  + "  len LF : " + str( len(LF[0])))
print(" mean LF[1]: " + str( np.mean( LF[1]) ).rjust(20) + "  var LF[1]: " +  str( np.var(LF[1])).rjust(20)  + "  len LF : " + str( len(LF[1])))
print(" mean LF[2]: " + str( np.mean( LF[2]) ).rjust(20) + "  var LF[2]: " +  str( np.var(LF[2])).rjust(20)  + "  len LF : " + str( len(LF[2])))
print(" mean LF[3]: " + str( np.mean( LF[3]) ).rjust(20) + "  var LF[3]: " +  str( np.var(LF[3])).rjust(20)  + "  len LF : " + str( len(LF[3])))
print("")



print("------- HF: difference function Y_l:")
print(get_Y_l_array(HF))
print("")

print("------- LF: difference function Y_l:")
print( get_Y_l_array(LF))
print("")

print("------- HF Var(Y_l_HF) :")
print( get_variances_for_diffs_Y_l( HF ) )
print("")
print("------- HF Var(Y_l_LF) :")
print( get_variances_for_diffs_Y_l( LF ) )
print("")


print("...........MLMC HF #samples: ...............................................................................................")

print( get_N_l_MLMC1(EPSILON_MLMF, 0,get_variances_for_diffs_Y_l(HF), COST_l_HF) )
print( get_N_l_MLMC1(EPSILON_MLMF, 1,get_variances_for_diffs_Y_l(HF), COST_l_HF) )
print( get_N_l_MLMC1(EPSILON_MLMF, 2,get_variances_for_diffs_Y_l(HF), COST_l_HF) )
print( get_N_l_MLMC1(EPSILON_MLMF, 3,get_variances_for_diffs_Y_l(HF), COST_l_HF) )

print("...........MLMF HF/LF: ...............................................................................................")


print("")
print("Correlations:")
print( get_RHO_l_s(HF,LF) )
print("")


print("LAMDBA: ")
print( lambda_aval_arr )
print("")
print("R_STAR: ")
print( r_values_list ) 
print("")


print("MLMF  HF #samples :")

print( get_N_l_MLMF(EPSILON_MLMF, 0, get_variances_for_diffs_Y_l(HF), CORRELATIONS_RHO_l, COST_l_HF) )
print( get_N_l_MLMF(EPSILON_MLMF, 1, get_variances_for_diffs_Y_l(HF), CORRELATIONS_RHO_l, COST_l_HF) )
print( get_N_l_MLMF(EPSILON_MLMF, 2, get_variances_for_diffs_Y_l(HF), CORRELATIONS_RHO_l, COST_l_HF) )
print( get_N_l_MLMF(EPSILON_MLMF, 3, get_variances_for_diffs_Y_l(HF), CORRELATIONS_RHO_l, COST_l_HF) )
print("")
print("MLMF LF #samples :")

for k in range(len(r_values_list)):
	print( "N_l_LF: " + str(n_samples_mlmf_var[k]*(1+r_values_list[k]) ) )

print("...............................................................................................")


print("FINAL RESULT  for QoI: "+ val_QoI +"   with epsilon: "+ str(math.sqrt(2/EPSILON_MLMF)) )
print("MLMC:           "+ str(result_mlmc_var) + "  MLMF:   " +str( result_mlmf_var ) + "  Benchmark: " + str( np.var(HF[3])) )
print("comp_cost_MLMC: "+ str(comp_cost_MLMC)  + "  comp_cost_MLMF: "+ str(comp_cost_MLMF) + " compared COST: "  + str( COST_l_HF[3]*len(HF[3]) ) )
