#!/bin/bash
# $Id$
# 
echo "Start cleanup"
make veryclean UNAME=Linux-gcc-8.3.1 BUILD_DIR=Linux-gcc-8.3.1 
echo "All done"
