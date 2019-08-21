# Compiling Command Line version of KMC1D on CUDA-capable GPU Linux sever

After downloading the source code of **KMC1D**, you need to ensure that you 1)have a CUDA-Capable **GPU Environment**; 2) install **NVIDIA** and **CUDA** related Toolkit or Packages; and 3) **Compile** the Source Code into executable binary command. 

## CUDA-capable GPU  Environment
**CUDA** is a **parallel** computing platform and programming model invented by **NVIDIA.** Essentially, NVIDIA is committed to supporting CUDA as hardware changes. Hardware is projected to change radically in the future. However, program algorithm, architecture and source code can remain largely unchanged. To use CUDA on your system, you need a CUDA-capable GPU sever, supported version of **Linux** with a **gcc** compiler and **nvcc** compiler, and NVIDIA CUDA related Toolkit. The details for CUDA installation guide can be found at [here](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html#pre-installation-actions) . Briefly, the main command lines for CUDA Installation are the followings:
### Verify CUDA-Capable GPU in your Linux sever
``` 
lspci | grep -i nvidia 
```
If your sever is equipped with GPU hardware, you will get the followings: 
```
00:07.0 3D controller: Corporation GV100GL [Tesla V100 PCIe 32GB] (rev a1)
00:08.0 3D controller: NVIDIA Corporation GV100GL [Tesla V100 PCIe 32GB] (rev a1)
```
### Install NVIDIA and CUDA related Toolkits or Packages
To **Ubuntu** or **Debian** like Linux, you only need to use "apt install" to install CUDA NVIDIA Toolkit
```
sudo apt install nvidia-cuda-toolkit 
``` 
To **Centos** or **Redhat** like Linux, you need to use "yum install" to install the related packages and set the CUDA Environment variables 
```
yum install epel-release*
yum install *cuda*
yum install *nvdia*
```
and 
### Set Environment Variables for CUDA
```
echo 'export PATH=/usr/local/cuda-9.0/bin${PATH:+:${PATH}}' >> ~/.bash_profile
echo 'export LD_LIBRARY_PATH=/usr/local/cuda 9.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}' >> ~/.bash_profile
```
After set the Environment Variables, check your .bash_profile file
```
less ~/.bash_profile
``` 
you should get the following display  
```
# .bash_profile
# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi
# User specific environment and startup programs
PATH=$PATH:$HOME/.local/bin:$HOME/bin
export PATH
export PATH=/usr/local/cuda-9.0/bin${PATH:+:${PATH}}
export LD_LIBRARY_PATH=/usr/local/cuda-9.0/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
```
## Compile the CUDA Source Code -KMC1D.cu
After download the KMC1D source code from [here](http://bioinfo.noble.org/KMC1D/download)
Use the following command lines to compile 
```
cd KMC1D
cd Source_Code
make clean
make
```
if compile successfully, you should get the following result
```
ls
KMC1D KMC1D.cu KMC1D.o Makefile
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTA5MDUxMjMyNSwtMTQ4NzMxNDg4OCwtMT
M0NzQ5MzcsMTUxODExMDMxNSwxNzM5NjQ1Nzc2LDEwNjMwMDA3
NjcsMTkzODM1OTcyOCwtMTM5MjUwNDI0NCwtMTI4NjMyODM3Ni
wtMTEwMDgxMTQ4XX0=
-->