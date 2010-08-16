#!/bin/bash
# All of the files for this lab can be obtained by direct download from
# PDC's Web site.  To simplify getting these files, copy this file 
# to your computer and execute it by typing:
#
# bash get_all.bash
# 
# It will download all the files from the Web for you.
# 
# The Mac does not have a built-in wget, so use curl instead by commenting
# in that line below and commenting out the wget one.

downloader="wget -c"
# downloader="curl -O"

base=http://www.pdc.kth.se/publications/talks/mpi/basics-lab

for i in hello.c \
    hello.ex1.c \
    hello.ex1.f \
    hello.ex1.f77 \
    hello.ex2.c \
    hello.ex2.f \
    hello.ex2.f77 \
    hello.f \
    hello.f77 \
    karp.c \
    karp.f \
    karp.f77 \
    karp.soln.c \
    karp.soln.f \
    karp.soln.f77 \
    data_values
do ${downloader} ${base}/$i
done
