# SemesterThesis
Semester Thesis: On Importance Sampling for Surrogate Model Construction

##########################################################################################################################
# FILES NEEDED:

In the directory are a number of json files, incl. a directory LF-200/ with json files.
These are the data input files from the Paul Scherrer Institut.
These files should not be changed.

The main driver file is extract_result.py


##########################################################################################################################
# HOW TO RUN?

There are 2 ways to run the calculation

1) Without Parameter

python3 extract_result.py  


2) With precision parameter epsilon, -e <VALUE> , where <VALUE> MUST be a number, for example like

python3 extract_result.py  -e 0.10697266120665758


The parameter -e <VALUE> corresponds to the epslion value in the script in the scetion 3.4 and 5.5.
	

##########################################################################################################################
# DISPLAYED RESULTS:

The displayed results will look like this

HF array 
 mean  HF[0]:    2.251186414356308  var HF[0]: 0.016671419475984953   len HF : 100
 mean  HF[1]:    2.166283214648033  var HF[1]: 0.006227494136198145   len HF : 100
 mean  HF[2]:   2.1774147040405536  var HF[2]: 0.014265634777419035   len HF : 100
 mean  HF[3]:    2.240524446535444  var HF[3]: 0.008505873688985748   len HF : 100

LF array 
 mean LF[0]:    3.003251973872809  var LF[0]:  0.04643285870897073  len LF : 496
 mean LF[1]:    2.350671683884118  var LF[1]: 0.008460076233196095  len LF : 500
 mean LF[2]:   2.2180571852876274  var LF[2]: 0.007651788460016149  len LF : 497
 mean LF[3]:   2.1511737378364852  var LF[3]: 0.007897851380519757  len LF : 495

------- HF: difference function Y_l:
[2.251186414356308, -0.08490319970827498, 0.01113148939252051, 0.06310974249489032]

------- LF: difference function Y_l:
[3.003251973872809, -0.6525802899886912, -0.13261449859649055, -0.06688344745114216]

------- HF Var(Y_l_HF) :
[0.016671419475984953, 0.01754582690973573, 0.018348491049883214, 0.022245669904917324]

------- HF Var(Y_l_LF) :
[0.04643285870897073, 0.05363998790229491, 0.01627120323636623, 0.01658431126114655]

...........MLMC HF #samples: ...............................................................................................
111.98344170924153
40.61715014050927
14.685127166376454
5.716830097160552
...........MLMF HF/LF: ...............................................................................................

Correlations:
[-0.17444933822313347, -0.33012097365557935, -0.413881067667887, 0.4390950312251155]

LAMDBA: 
[1.0015756902349138, 0.9490859722415099, 0.8989091733774546, 0.8807064399063641]

R_STAR: 
[-0.049227610131445276, 0.8768328672252077, 1.4399026017417897, 1.6228012749141603]

MLMF  HF #samples :
100.00000000000001
36.2706749502925
13.11365943235209
5.105067329511062

MLMF LF #samples :
N_l_LF: 95.07723898685549
N_l_LF: 69.44281608733269
N_l_LF: 34.15863642438506
N_l_LF: 15.736807649484962
...............................................................................................
FINAL RESULT  for QoI: HALO_X   with epsilon: 0.10697266120665758
MLMC:           0.005721553113009385  MLMF:   0.005060478069664908  Benchmark: 0.008505873688985748
comp_cost_MLMC: 123949.50766612907  comp_cost_MLMF: 136577.21801538247 compared COST: 1474560.0
