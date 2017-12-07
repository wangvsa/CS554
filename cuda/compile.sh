#!/bin/bash

nvcc test.cu ../matrix_io.cc ../mmio.cc ../util.cc --std=c++11
