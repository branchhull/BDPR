# BDPR

# Blind Deconvolutional Phase Retrieval (BDPR)
This page contains the different versions of BDPR. The first version was presented at NeurIPS 2018 by Ahmed, Aghasi and Hand:

* [**Link to NeurIPS paper**](http://papers.nips.cc/paper/8207-blind-deconvolutional-phase-retrieval-via-convex-programming)

The Matlab code for the first version is available in [*BDPR BDPR (First Version 2018)*]

# Faster Matlab Code (Latest Version) and Installation Guide
Later in 2019, the authors presented a new algorithm which uses a low rank factorization and speeds up the solution to the BDPR semidefinite program. The new paper which includes the new algorithm and the stability analysis is available at:

* [**Link to 2019 paper**](https://arxiv.org/)

The latest BDPR code uses the [**minFunc**](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) toolbox for intermediate nonlinear solves, which needs to be downloaded separately. To install the BDPR latest version please follow the easy steps below:

* download and unzip [*BDPR (Latest Version 2019)*]
* download [*minFunc.zip*] from [**here**](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html) and save it in the [*BDPR (Latest Version 2019)*] folder
* open [*BDPR (Latest Version 2019)*] in Matlab and run [*Setup_BDPR_minFunc.m*] (please make sure that mex is working properly in your Matlab for flawless compilation of the minFunc .C files)
* now you are ready to run the BDPR demo file

The file [*Setup_BDPR_minFunc.m*] contains a sequence of simple commands to unzip the downloaded folder, compile it and add it to the path. The user can easily apply those steps manually




