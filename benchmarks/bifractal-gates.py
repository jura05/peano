#!/usr/bin/env python3
"""
Benchmark: find possible gates for 2d bifractals of genus 4.

Time: 2m45s
Result: 9 good gates
"""

import sys
sys.path.append('.')

from peano.curves import gen_possible_gates


if __name__ == "__main__":
    list(gen_possible_gates(dim=2, div=2, pattern_count=2))
