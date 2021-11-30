# AbCD

## Ab-initio Catalyst Dynamic

including feature
-----------------
1. Mean-field Microkinetic Modeling
2. Parameter Estimation
3. Bayesian Inference

required packages
-----------------
1. numpy
2. scipy
3. casadi_2.4.2 (https://github.com/casadi/casadi/releases?after=3.0.0-rc3)


## Installation
1. install python=2.7 and required packages (using conda to create myenv environment)
```
conda create -n myenv python=2.7
conda install -n myenv numpy scipy
```
2. install casadi_2.4.2 from (https://github.com/casadi/casadi/releases/tag/2.4.2)
    add to myenv conda environment
3. install AbCD
```
git clone https://github.com/thj2009/AbCD.git
conda develop AbCD
```

4. test the installation
```
python testing/test_cstr.py
```
