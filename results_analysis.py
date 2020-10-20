import pickle
from statistics import mean, pstdev
from math import isclose

from dwave.embedding import unembed_sampleset

def HammingDistance(vec1, vec2):
    return float(sum(v1 ^ v2 for v1, v2 in zip(vec1, vec2)))

# General parameters
gsFractions = {}
hDistances = {}
numSamples = 1000
nurses = 3

for days in range(5, 15):

    # Expected ground state energy
    gs = 0.3 * (days % nurses) 

    # Output data for analysis
    filename = 'reverse_results_chimera_N%d_D%d_s%d.p' % (nurses, days, numSamples)
    saveData = pickle.load(open(filename, "rb"))

    results = saveData['results']
    embedding = saveData['embedding']
    bqm = saveData['bqm']

    # If unembedded results are not available for some reason, catch exception
    try:
        samples = saveData['samples']
    except KeyError:
        samples = unembed_sampleset(results, embedding, bqm, chain_break_fraction=True)

    # print some basic information to see what we are doing if needed
    #print(results.info)
    #print(embedding)
    #print(samples.record[0:10])
    #print(samples.first)

    # Lambda function to calculate index(n, d) as a function of nurse and day.
    index = lambda n, d: n * days + d

    # Print first schedule obtained (idealy, a ground state)
    print('Solution with energy = {}.'.format(samples.first.energy) )

    header = " nurse | "
    for day in range(days):
        header += "%02d | " % day
    print(header)

    for nurse in range(nurses):
        s = "   %d   | " % (nurse + 1)
        for day in range(days):
            s += ("   | ", "XX | ")[samples.first.sample[index(nurse,day)]]
        print(s)

    ## Check Hamming distances between consecutive states
    distances = []
    for i in range(len(samples.record.sample) - 1):
        h = HammingDistance(samples.record.sample[i], samples.record.sample[i + 1])
        distances.append(h)

    mean_h = mean(distances)
    std_h = pstdev(distances)

    hDistances[nurses, days] = (mean_h, std_h)

    # Fraction of solutions in ground state
    numGs = sum(r.num_occurrences for r in samples.record if isclose(r.energy, gs, abs_tol=1e-3) )
    p = numGs / numSamples
    gsFractions
[nurses, days] = p

print(gsFractions)
print(hDistances)