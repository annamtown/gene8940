#!/bin/bash

cd /escratch4/s_11/s_11_Aug_17/fasttree

#path for fasttree
export PATH=/usr/local/fasttree/2.1.4/:$PATH

FastTree -nt -gtr < /escratch4/s_11/s_11_Aug_17/project2/allconsensus.phy > fasttree_file
