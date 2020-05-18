#!/bin/env python3

from collections import Counter
from math import log

str = "AAAATTTGGGG"
l = len(str)
counts = Counter(str)
print(counts)

parts = [(i/l)*(log((i/l),2)) if i  != 0 else 0 for i in counts.values()]

print(parts)
Shannon = sum(parts)
print(-Shannon)
