import numpy as np
from Extractables.InductionModel import LayeredSystem

test = LayeredSystem(json=r"./Working/test.json")

print(test._structure)
print(test._name)
print(test._maxlayers)