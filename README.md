# Power Flow

This repository contains the implementation of power flow analysis for various systems, including 2-bus, 33-bus, 69-bus, and 118-bus networks. It also provides multiple custom power flow algorithms used in the accompanying research paper.

## Folder Structure

- **data**: This folder includes the test case data for different systems used in the paper. The available systems are:
  - 2-bus system
  - 33-bus system
  - 69-bus system
  - 118-bus system

- **lib**: This folder contains the necessary functions from the MATPOWER library, as well as custom power flow algorithms implemented in the paper. The algorithms provided are:
  - **NR**
  - **CJ**
  - **CJ3-Vk**
  - **CJ3-Vtemp**
  - **CJ2-1**
  - **CJ2-2**
  - **FDXB**

## Usage

To run the power flow analysis, follow these steps:

1. Open MATLAB and navigate to the `power-flow` directory.

2. Set the power flow algorithm and run the power flow analysis using the following command:
   ```matlab
   mpopt = mpoption('pf.alg','algorithm'); 
   results = runpf('casedata', mpopt);
3. Replace `'algorithm'` with the desired algorithm name (e.g., `'NB'`, `'CJ'`, etc.) and `'casedata'` with the desired case file (e.g., `'case2'`, `'case33'`, `'case69'`, or `'case118'`).

    **Example**:
    ```matlab
    mpopt = mpoption('pf.alg','NR'); 
    results = runpf('case2', mpopt);