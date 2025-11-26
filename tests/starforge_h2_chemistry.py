from jaco.models.starforge import h2_chemistry

for p in h2_chemistry.h2_chemistry.subprocesses:
    print(p.equation, p.rate, p.bibliography)
