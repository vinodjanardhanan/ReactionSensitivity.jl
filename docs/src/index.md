```@meta
CurrentModule = ReactionSensitivity
```

# ReactionSensitivity
The sensitivity analysis evaluates the global sensitivity using Sobol sampling.  It makes use of GlobalSensitivity.jl package for the evaluation of sensitivity indices. It reports the Sobol indices for all species mentioned in the input file for the micro-kinetic model parameters. The function evaluation may be performed using a batch reactor model, stirred reactor model or plug flow reactor model.  You may please refer to the documentation for the respective reactor models to understand the governing equations that are solved. 

In general, only the forward reaction rate parameters can be specified independently for a micro-kinetic model where the reactions are expressed as reversible pairs. The reverse reaction rate parameters are then calculated from the equilibrium constant. However, this requires thermodynamic data for the calculation of equilibrium constants. 
In many surface reaction mechanisms, the reversible reactions are explicitly expressed as two reactions; one represents the forward direction, and the other represents the reverse direction. This difficulty is mainly due to the non-availability of the NASA polynomial coefficients for the surface adsorbed species to calculate equilibrium constants. However, for such a reaction mechanism to be thermodynamically consistent, one has to ensure the following constraint on all forward-reverse reaction pairs.

```math
\frac{A_fT^\beta}{A_r T^\beta} = \exp(\frac{\Delta S}{R})
```
In such a case, all parameters correlate to at least one parameter. Therefore altering one parameter will alter another parameter constraint to the above correlation. While developing a surface reaction mechanism, one may assume that all parameters are uncorrelated. For instance, if the reaction mechanism contains 20 reversible reactions expressed as 40 irreversible pairs, one may want to perform a full sensitivity analysis assuming that all 40 parameters are uncorrelated. At some stage during the mechanism development, if you make the mechanism thermodynamically consistent, you may prefer to perform the analysis using only 20 independent parameters. In that case, altering any of the 20 parameters will also alter the parameter of the corresponding reaction pair. One can also think of secondary interactions, which are not considered here. So, the analysis assumes that the parameters are independent or uncorrelated in either case. This is important for performing Sobol decomposition or ANOVA decomposition. 


Documentation for [ReactionSensitivity](https://github.com/vinodjanardhanan/ReactionSensitivity.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("ReactionSensitivity")
```


## General interfaces
```@index
```

```@autodocs
Modules = [ReactionSensitivity]
```

## Theory
The Sobol decomposition of a function $y=f(x)$ is given as

```math
y= f(x) = f_0 + \sum_{i=1}^n f_i(x_i) + \sum_{i=1}^n \sum_{j=i+1}^n f_{ij}(x_i,x_j) + \ldots + f_{12\cdots n}(x_1, x_2,\ldots,x_n)
```
where
```math
f_0 = E(y)
```

```math
f_i(x_i) = E(y\vert x_i) - E(y)
```

```math
f_{ij}(x_i,x_j) = E(y\vert x_i,x_j) - E(y\vert x_i) - E(y\vert x_j) - E(y) = E(y\vert x_i,x_j) - f_i - f_j - f_0
```
Taking the variance of function $f(x)$

```math
V(y) = V[f_0] +V\left[ \sum_{i=1}^n f_i(x_i)\right] + V\left[\sum_{i=1}^n \sum_{j=i+1}^n f_{ij}(x_i,x_j)\right] + \ldots +V[ f_{12\cdots n}(x_1, x_2,\ldots,x_n)]
```
i.e.,

```math
V(y) = 0 +  \sum_{i=1}^n V[f_i(x_i)] + \sum_{i=1}^n \sum_{j=i+1}^n V[f_{ij}(x_i,x_j)] + \ldots + V[f_{12\cdots n}(x_1, x_2,\ldots,x_n)]
```

```math
V = \sum_{i=1}^n V_i + \sum_{i=1}^n \sum_{j=i+1}^n V_{ij}+ \ldots + V_{12\ldots n}
```

The global sensitivity index is then defined as

```math
S_i =\frac{ V_i}{V}
```
which is also known as the main sensitivity index of the first order Sobol index. The total sensitivity index is defined as

```math
S_i^T = \frac{V_i^T}{V}
```
where
```math
V_i^T= V_i + \sum_{j} V_{ij} + \sum_{jk} V_{ijk} + \ldots
```
The first order Sobol index gives the variance of the model output due to the corresponding parameter input. The total sensitivity index also considers the interaction among the parameters. Therefore, the difference between the total and main sensitivity indices shows the contribution due to interaction. 

## Running the code
On the Julia REPL 
```julia
julia>using ReactionSensitivity
julia>rxn_gsa("sensitivity.xml","../lib/")
```

## Input file
The above method takes two input arguments. The first argument points to an xml input files that defines the species and parameters for which the sensitivity analysis needs to be performed. The second argument points to the directory in which the "therm.dat" and the mechanism input files are present. 

The structure of the XML input file is shown below.
```
?xml version="1.0" encoding="ISO-8859-1"?>
<sensitivity>
    <gsa_species>CH4 H2O H2 CO CO2 O2 N2</gsa_species>
    <gsa_lb>30572</gsa_lb>
    <gsa_ub>10</gsa_ub>
    <gsa_N>100</gsa_N>
    <gsa_model>cstr.xml</gsa_model>
    <p_smech>1:6=7:12,15=16,16:20</p_smech>
</sensitivity>
```
The meaning of different tags is listed below.

- <sensitivity>: Root xml tag for sensitivity analysis
- <gsa\_species>: The species that needs to be monitored for the parameter sensitivity analysis. The code will produce the output of parameter sensitivity of only those species which are listed in this element
- <gsa\_lb>: lower percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the lower bound considered for design matrix generation is 90 with a gsa_lb = 10
- <gsa\_ub>: upper percentage bound for the parameter value. For instance, if the value of the actual parameter is 100, then the upper bound considered for design matrix generation is 110 with gsa_ub = 10
- <gsa\_N>: the size of the sample size. With a sample size of 100, 100 points will be sampled between the lower and upper bounds for each parameter. The higher the sample size higher will be calculation time. (suggested : 2^p, where p is the number of parameters)
- <gsa\_model>: the reactor model that will be considered for integration. The content of the tag is the name of the XML input file of the reactor model. The file must be present in the working directory. You may copy the input file for Batch reactor, StirredReactor or PlugFlowReactor. 
- <p\_smech>: parameters of the surface reaction mechanism that needs to be considered for sensitivity analysis. There are different options for specifying the parameter list from a mechanism
    -   All parameters: Use ":" as the value of <p\_smech> element to consider all pre-exponential factors and sticking coefficient as sensitivity parameters. In such specifications, the reverse reaction rates will not be adjusted to maintain the initial ratio. (e.g. <p\_smech>:</p_smech>)
    -   Reversible reaction: If the reverse reaction parameters need to be adjusted according to the initial ratio, then a reverse reaction must be specified for each forward reaction. For instance, if you want to consider the pre-exponential factor of reaction 15 as a sensitivity parameter, during sensitivity analysis, this parameter will be sampled between its upper and lower limit. This also changes the ratio between the forward and reverse reaction rate parameters as specified in the original mechanism. By specifying a reverse reaction, the ratio can be maintained a constant. For instance, in the above example, reaction 16 is the reverse of reaction 15. For every sample values of the pre-exponential factor of reaction 15, the pre-exponential factor of reaction 16 will be adjusted to maintain the initial ratio between these two in the mechanism input file. Only the parameter of reaction to the left of "=" will be considered as sensitivity parameter, and the pre-exponential factor to the right of "=" will be adjusted to maintain the initial ratio. 
    -   Contiguous parameters: Use ":" to specify contiguous reactions. For instance, 1:6=7:12 means pre-exponential factors of reactions from 1 to 6 will be considered as sensitivity parameters, and reactions from 7 to 12 are reverse pairs of reactions from 1 to 6. i.e., reaction 7 is the reverse pair of reaction 1, reaction 8 is the reverse pair of reaction 2 and so on. Comma-separated list of parameters can be provided in <p\_smech>
    -   Unconstrained reverse reactions: If the user does not want to constrain the reverse reaction parameters according to the initial ratio, then omit "=". For instance, <p\_smech>16:20</p_smech> means the pre-exponential factors of reactions 16 to 20 will be considered as sensitivity parameters. However, their corresponding reverse reaction parameters will not be adjusted to maintain the initial ratio.

## Input file download
The xml input file and the *lib* directory containig other required input files may be downloaded from [here](https://github.com/vinodjanardhanan/ReactionSensitivity.jl/tree/main/test).



## Output
### Sobol first order indices
```
 parameter             CH4             H2O              H2              CO             CO2              O2              N2
       k13      +3.4690e-08     -1.7068e-08     -2.8338e-09     -1.9325e-09     -1.7134e-08     -1.3293e-08     +1.1196e-05
       k15      -6.5332e-03     +3.3629e-03     +1.8052e-03     +1.8235e-03     +3.4042e-03     -5.7777e-04     -1.3007e-01
       k17      -3.9335e-04     +8.3314e-05     +1.7312e-05     +1.7862e-05     +8.3680e-05     -9.3703e-05     -6.8455e-03
       k19      -6.5291e-08     -6.0111e-08     +6.8006e-09     +8.0495e-09     -6.5388e-08     -6.8575e-09     -1.7432e-05
       k21      +7.3879e-04     -6.6383e-04     -4.2291e-04     -4.2971e-04     -6.6497e-04     -7.9864e-05     +1.7248e-02
       k23      -4.8677e-08     -1.6787e-08     -4.9927e-09     -3.1755e-09     -9.9494e-09     -3.2560e-08     -4.4434e-06
       k25      -1.7704e-08     +1.0263e-08     -4.4179e-09     -3.0041e-09     +1.8139e-08     -7.4624e-09     -1.4354e-06
       k27      +1.1487e-03     +2.4962e-04     +3.6667e-03     +3.5571e-03     +6.5261e-04     -9.5928e-07     +1.2910e-02
       k29      -2.1521e-03     +1.8689e-03     +1.4845e-03     +1.4913e-03     +1.8540e-03     +2.8824e-04     -5.4827e-02
       k31      -2.7705e-04     +1.6704e-04     +1.3851e-04     +1.3912e-04     +1.6336e-04     -5.0479e-07     -6.4453e-03
       k33      -1.2737e-07     +5.8412e-09     -1.3247e-07     -1.2771e-07     +3.3854e-08     -6.1005e-09     -3.0165e-05
       k35      +9.6615e-01     +9.8403e-01     +9.6273e-01     +9.6294e-01     +9.8475e-01     +9.7775e-01     +8.4048e-01
       k37      -1.4904e-03     +1.6085e-03     +9.7492e-04     +9.8230e-04     +1.6416e-03     +3.5735e-04     -3.8559e-02
       k39      +2.7518e-04     -1.6529e-04     -1.2171e-04     -1.2261e-04     -1.6221e-04     +1.5000e-07     +6.2832e-03
       k41      -2.9709e-08     -7.7295e-09     -1.7599e-09     -6.4605e-10     -2.6278e-09     -1.6481e-08     -2.9450e-06

```
### Sobol total indices
```
 parameter             CH4             H2O              H2              CO             CO2              O2              N2
       k13      -2.1533e-08     -2.7865e-08     -4.4028e-09     -3.2608e-09     -2.7959e-08     +1.6423e-09     -1.0754e-05
       k15      +1.4674e-02     +1.9171e-03     +8.6542e-03     +8.6329e-03     +1.3887e-03     +6.9923e-03     +1.4053e-01
       k17      +5.7453e-04     +2.1747e-06     +2.5010e-04     +2.4889e-04     -1.2390e-05     +2.1545e-04     +7.1129e-03
       k19      +6.4776e-08     +5.6530e-08     -1.5435e-08     -1.6307e-08     +6.1020e-08     +5.7361e-09     +1.7467e-05
       k21      +1.1296e-04     +1.1520e-03     +1.5825e-03     +1.5845e-03     +1.1032e-03     +7.1078e-04     -1.6090e-02
       k23      +5.8070e-08     -1.4402e-08     -1.8386e-09     -1.9859e-09     -2.0878e-08     +2.4275e-08     +4.7889e-06
       k25      +3.5029e-08     -1.8926e-08     +1.1639e-10     -8.8117e-12     -2.4413e-08     +1.5571e-08     +1.7272e-06
       k27      +1.5784e-03     +1.1853e-03     +7.5254e-03     +7.3373e-03     +1.8150e-03     +4.3244e-04     -1.8205e-03
       k29      +1.2231e-02     +4.8143e-03     +1.1300e-02     +1.1290e-02     +4.2413e-03     +7.7509e-03     +6.7611e-02
       k31      +3.7639e-04     -1.1361e-04     +2.8578e-07     -5.9457e-07     -1.1713e-04     +7.1714e-05     +6.5840e-03
       k33      +9.0423e-08     +5.6970e-08     +3.2829e-08     +2.9780e-08     +4.8873e-08     +2.4054e-08     +3.0129e-05
       k35      +1.0056e+00     +9.9598e-01     +9.9511e-01     +9.9520e-01     +9.9543e-01     +1.0005e+00     +1.1175e+00
       k37      +2.3898e-02     +1.4246e-02     +2.6681e-02     +2.6668e-02     +1.3084e-02     +1.8109e-02     +6.6213e-02
       k39      -2.7029e-04     +1.6445e-04     +1.3241e-04     +1.3296e-04     +1.6092e-04     +1.1004e-06     -6.2724e-03
       k41      +4.1251e-08     -7.1443e-09     -2.5744e-09     -2.5212e-09     -1.1046e-08     +1.7981e-08     +3.1932e-06

```