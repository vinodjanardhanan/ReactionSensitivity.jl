var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ReactionSensitivity","category":"page"},{"location":"#ReactionSensitivity","page":"Home","title":"ReactionSensitivity","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The sensitivity analysis evaluates the global sensitivity using Sobol sampling.  It makes use of GlobalSensitivity.jl package for the evaluation of sensitivity indices. It reports the Sobol indices for all species mentioned in the input file for the micro-kinetic model parameters. The function evaluation may be performed using a batch reactor model, stirred reactor model or plug flow reactor model.  You may please refer to the documentation for the respective reactor models to understand the governing equations that are solved. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"In general, only the forward reaction rate parameters can be specified independently for a micro-kinetic model where the reactions are expressed as reversible pairs. The reverse reaction rate parameters are then calculated from the equilibrium constant. However, this requires thermodynamic data for the calculation of equilibrium constants.  In many surface reaction mechanisms, the reversible reactions are explicitly expressed as two reactions; one represents the forward direction, and the other represents the reverse direction. This difficulty is mainly due to the non-availability of the NASA polynomial coefficients for the surface adsorbed species to calculate equilibrium constants. However, for such a reaction mechanism to be thermodynamically consistent, one has to ensure the following constraint on all forward-reverse reaction pairs.","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracA_fT^betaA_r T^beta = exp(fracDelta SR)","category":"page"},{"location":"","page":"Home","title":"Home","text":"In such a case, all parameters correlate to at least one parameter. Therefore altering one parameter will alter another parameter constraint to the above correlation. While developing a surface reaction mechanism, one may assume that all parameters are uncorrelated. For instance, if the reaction mechanism contains 20 reversible reactions expressed as 40 irreversible pairs, one may want to perform a full sensitivity analysis assuming that all 40 parameters are uncorrelated. At some stage during the mechanism development, if you make the mechanism thermodynamically consistent, you may prefer to perform the analysis using only 20 independent parameters. In that case, altering any of the 20 parameters will also alter the parameter of the corresponding reaction pair. One can also think of secondary interactions, which are not considered here. So, the analysis assumes that the parameters are independent or uncorrelated in either case. This is important for performing Sobol decomposition or ANOVA decomposition. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for ReactionSensitivity.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use the following commands in the julia REPL","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"ReactionSensitivity\")","category":"page"},{"location":"#General-interfaces","page":"Home","title":"General interfaces","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ReactionSensitivity]","category":"page"},{"location":"#ReactionSensitivity.get_rxn_ids-Tuple{LightXML.XMLElement, AbstractString}","page":"Home","title":"ReactionSensitivity.get_rxn_ids","text":"A function to extract the parameters for sensitivity analysis from the xml input      The function returns two integer arrays. forwardrxnids contains the ids of the      reactions whoes parameters will be considered as sensitivity parameters.     reverserxnids contains the ids of the reactions aginst with the forward reaction      parameters are correlated. \n\n\n\n\n\n","category":"method"},{"location":"#ReactionSensitivity.gsa_problem-Tuple{Any, Any, Any, ReactionSensitivity.Gsa_parameters, ReactionSensitivity.Solver}","page":"Home","title":"ReactionSensitivity.gsa_problem","text":"A common function for all reactor models to define the sensitivity problem \n\nUsage\n\ngsa_problem(params, prob, gsa_p)\n\nparams : a tuple of parameters returned by the respective reactor model invoked \nprob : the problem definition for the reactor model \ngsap : struct of the type gsaparameters\n\n\n\n\n\n","category":"method"},{"location":"#ReactionSensitivity.gsa_smech_parameter_map!-Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Any}","page":"Home","title":"ReactionSensitivity.gsa_smech_parameter_map!","text":"This function returns the sticking coefficient or the pre-exponent factor      corresponding to the reaction ids present in gsarxnids \n\nUsage\n\ngsa_smech_parameter_map!(gsa_smech_params, gsa_srxn_ids, gsa_srxn_constraint_ids, md )\n\ngsasmechparams : vector containing the parameters extracted from the inoput mechanism (output)\ngsasrxnids : reaction ids whoes parameters needs to be extracted from the mechanism (input)\ngsasrxnconstraint_ids : id of the reactions which needs to be considered for reverse reaction \nmd: MechanismDefinition\n\n\n\n\n\n","category":"method"},{"location":"#ReactionSensitivity.initial_parameter_ratio!-Tuple{Vector{Float64}, Vector{Int64}, Vector{Int64}, Any}","page":"Home","title":"ReactionSensitivity.initial_parameter_ratio!","text":"A function to calculate the intial ratio of pre-exponential factors\n\n\n\n\n\n","category":"method"},{"location":"#Theory","page":"Home","title":"Theory","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Sobol decomposition of a function y=f(x) is given as","category":"page"},{"location":"","page":"Home","title":"Home","text":"y= f(x) = f_0 + sum_i=1^n f_i(x_i) + sum_i=1^n sum_j=i+1^n f_ij(x_ix_j) + ldots + f_12cdots n(x_1 x_2ldotsx_n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"where","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_0 = E(y)","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_i(x_i) = E(yvert x_i) - E(y)","category":"page"},{"location":"","page":"Home","title":"Home","text":"f_ij(x_ix_j) = E(yvert x_ix_j) - E(yvert x_i) - E(yvert x_j) - E(y) = E(yvert x_ix_j) - f_i - f_j - f_0","category":"page"},{"location":"","page":"Home","title":"Home","text":"Taking the variance of function f(x)","category":"page"},{"location":"","page":"Home","title":"Home","text":"V(y) = Vf_0 +Vleft sum_i=1^n f_i(x_i)right + Vleftsum_i=1^n sum_j=i+1^n f_ij(x_ix_j)right + ldots +V f_12cdots n(x_1 x_2ldotsx_n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"i.e.,","category":"page"},{"location":"","page":"Home","title":"Home","text":"V(y) = 0 +  sum_i=1^n Vf_i(x_i) + sum_i=1^n sum_j=i+1^n Vf_ij(x_ix_j) + ldots + Vf_12cdots n(x_1 x_2ldotsx_n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"V = sum_i=1^n V_i + sum_i=1^n sum_j=i+1^n V_ij+ ldots + V_12ldots n","category":"page"},{"location":"","page":"Home","title":"Home","text":"The global sensitivity index is then defined as","category":"page"},{"location":"","page":"Home","title":"Home","text":"S_i =frac V_iV","category":"page"},{"location":"","page":"Home","title":"Home","text":"which is also known as the main sensitivity index of the first order Sobol index. The total sensitivity index is defined as","category":"page"},{"location":"","page":"Home","title":"Home","text":"S_i^T = fracV_i^TV","category":"page"},{"location":"","page":"Home","title":"Home","text":"where","category":"page"},{"location":"","page":"Home","title":"Home","text":"V_i^T= V_i + sum_j V_ij + sum_jk V_ijk + ldots","category":"page"},{"location":"","page":"Home","title":"Home","text":"The first order Sobol index gives the variance of the model output due to the corresponding parameter input. The total sensitivity index also considers the interaction among the parameters. Therefore, the difference between the total and main sensitivity indices shows the contribution due to interaction. ","category":"page"},{"location":"#Running-the-code","page":"Home","title":"Running the code","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"On the Julia REPL ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using ReactionSensitivity\njulia>rxn_gsa(\"sensitivity.xml\",\"../lib/\")","category":"page"},{"location":"#Input-file","page":"Home","title":"Input file","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The above method takes two input arguments. The first argument points to an xml input files that defines the species and parameters for which the sensitivity analysis needs to be performed. The second argument points to the directory in which the \"therm.dat\" and the mechanism input files are present. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The structure of the XML input file is shown below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<sensitivity>\n    <gsa_species>CH4 H2O H2 CO CO2 O2 N2</gsa_species>\n    <gsa_lb>5</gsa_lb>\n    <gsa_ub>10</gsa_ub>\n    <gsa_N>100</gsa_N>\n    <gsa_model>cstr.xml</gsa_model>\n    <p_smech>1:6=7:12,15=16,16:20</p_smech>\n</sensitivity>","category":"page"},{"location":"","page":"Home","title":"Home","text":"The meaning of different tags is listed below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<sensitivity>: Root xml tag for sensitivity analysis\n<gsa_species>: The species that needs to be monitored for the parameter sensitivity analysis. The code will produce the output of parameter sensitivity of only those species which are listed in this element\n<gsa_lb>: lower percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the lower bound considered for design matrix generation is 90 with a gsa_lb = 10\n<gsa_ub>: upper percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the upper bound considered for design matrix generation is 110 with gsa_ub = 10\n<gsa_N>: the size of the sample size. With a sample size of 100, 100 points will be sampled between the lower and upper bounds for each parameter. The higher the sample size higher will be calculation time. \n<gsa_model>: the reactor model that will be considered for integration. The content of the tag is the name of the XML input file of the reactor model. The file must be present in the working directory. You may copy the input file for Batch reactor, StirredReactor or PlugFlowReactor. \n<p_smech>: parameters of the surface reaction mechanism that needs to be considered for sensitivity analysis. There are different options for specifying the parameter list from a mechanism\nAll parameters: Use \":\" as the value of <p_smech> element to consider all pre-exponential factors and sticking coefficient as sensitivity parameters. In such specifications, the reverse reaction rates will not be adjusted to maintain the initial ratio. (e.g. <p_smech>:</p_smech>)\nReversible reaction: If the reverse reaction parameters need to be adjusted according to the initial ratio, then a reverse reaction must be specified for each forward reaction. For instance, if you want to consider the pre-exponential factor of reaction 15 as a sensitivity parameter, during sensitivity analysis, this parameter will be sampled between its upper and lower limit. This also changes the ratio between the forward and reverse reaction rate parameters as specified in the original mechanism. By specifying a reverse reaction, the ratio can be maintained a constant. For instance, in the above example, reaction 16 is the reverse of reaction 15. For every sample values of the pre-exponential factor of reaction 15, the pre-exponential factor of reaction 16 will be adjusted to maintain the initial ratio between these two in the mechanism input file. Only the parameter of reaction to the left of \"=\" will be considered as sensitivity parameter, and the pre-exponential factor to the right of \"=\" will be adjusted to maintain the initial ratio. \nContiguous parameters: Use \":\" to specify contiguous reactions. For instance, 1:6=7:12 means pre-exponential factors of reactions from 1 to 6 will be considered as sensitivity parameters, and reactions from 7 to 12 are reverse pairs of reactions from 1 to 6. i.e., reaction 7 is the reverse pair of reaction 1, reaction 8 is the reverse pair of reaction 2 and so on. Comma-separated list of parameters can be provided in <p_smech>\nUnconstrained reverse reactions: If the user does not want to constrain the reverse reaction parameters according to the initial ratio, then omit \"=\". For instance, <p_smech>16:20</p_smech> means the pre-exponential factors of reactions 16 to 20 will be considered as sensitivity parameters. However, their corresponding reverse reaction parameters will not be adjusted to maintain the initial ratio.","category":"page"},{"location":"#Input-file-download","page":"Home","title":"Input file download","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The xml input file and the lib directory containig other required input files may be downloaded from here.","category":"page"},{"location":"#Output","page":"Home","title":"Output","text":"","category":"section"},{"location":"#Sobol-first-order-indices","page":"Home","title":"Sobol first order indices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":" parameter             CH4              H2              CO\n        k1      +2.3448e-06     -1.4716e-05     +1.0065e-04\n        k2      -5.4885e-07     -1.3963e-05     +3.7032e-04\n        k3      -1.7705e-06     -6.2732e-06     +2.0993e-04\n        k4      -2.9290e-06     +3.2338e-06     +5.0191e-04\n        k5      -9.9537e-07     +2.8037e-05     -2.3794e-04\n        k6      +2.0481e-06     -1.7069e-05     +4.4028e-04\n       k13      -1.4126e-06     -4.9288e-06     +1.2184e-04\n       k15      +2.4661e-01     +8.8478e-01     +2.8693e-01\n       k17      -5.3015e-06     +7.5004e-06     +3.0034e-04\n       k19      +3.0071e-01     +8.0448e-02     +3.9064e-01\n       k21      +4.5341e-06     -1.1998e-05     +7.5951e-05\n       k23      -8.9618e-07     -7.9828e-06     +1.0737e-04\n       k25      +3.5734e-01     +1.4675e-01     +3.5949e-02\n       k27      +6.8595e-02     +1.7802e-01     +5.4256e-01\n       k29      -2.3420e-04     +3.0240e-03     +9.6439e-03\n       k31      +1.0861e-03     -5.9457e-05     -2.6299e-04\n       k33      +2.7207e-04     +5.4915e-04     +1.7431e-03\n       k35      -1.0602e-04     +6.1263e-04     +2.7361e-03\n       k37      +3.2406e-06     -1.9474e-05     +3.9189e-04\n       k39      -4.8068e-06     +1.3804e-05     +3.0318e-04\n       k41      -6.1879e-06     +1.7323e-05     +1.8869e-04","category":"page"},{"location":"#Sobol-total-indices","page":"Home","title":"Sobol total indices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":" parameter             CH4              H2              CO\n        k1      -2.7004e-06     +1.4402e-05     -1.0065e-04\n        k2      +5.2612e-07     +1.3897e-05     -3.6943e-04\n        k3      +1.3854e-06     +5.9647e-06     -2.0954e-04\n        k4      +2.8990e-06     -3.3247e-06     -5.0210e-04\n        k5      +5.4416e-06     -2.9566e-05     +3.5093e-04\n        k6      -2.2553e-06     +1.6854e-05     -4.3998e-04\n       k13      +1.1156e-06     +4.7096e-06     -1.2221e-04\n       k15      +3.9272e-01     -1.2662e-01     -1.9446e-01\n       k17      +5.1908e-06     -7.6530e-06     -2.9940e-04\n       k19      +2.5135e-01     +2.9465e-01     +9.3012e-01\n       k21      -4.1005e-06     +1.1996e-05     -7.8270e-05\n       k23      +7.3471e-07     +7.8210e-06     -1.0714e-04\n       k25      +2.8921e-01     +6.1753e-01     +6.1578e-02\n       k27      +1.1322e-01     -5.4445e-02     -5.0018e-02\n       k29      +8.3402e-04     -2.6450e-03     -7.7674e-03\n       k31      +8.5053e-04     +1.4107e-03     +5.0973e-03\n       k33      +4.0887e-04     -7.9329e-05     +1.0736e-06\n       k35      +1.3536e-04     -5.8887e-04     -2.6520e-03\n       k37      -3.1684e-06     +1.9450e-05     -3.9074e-04\n       k39      +4.7581e-06     -1.3864e-05     -3.0281e-04\n       k41      +6.2943e-06     -1.7289e-05     -1.8872e-04","category":"page"}]
}