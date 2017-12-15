#!/bin/bash

nvcc main.cu esc.cu ../matrix_io.cc ../mmio.cc ../util.cc --std=c++11
