# Nurse shift replication by Cassio Amorim, CJS Inc.
# Licensed under Apache 2.0 License.
# Original paper by Ikeda, Nakamura, Humble. DOI: 10.1038/s41598-019-49172-3
# Original licensed under Creative Commons

## Target Hamiltonian:
## 目的ハミルトニアン：
## H(q) = \sum_n,n' \sum_d, d' J_i(n, d)j(n',d')q_iq_j
##       + λ \sum_d [\sum_n E(n)*q_i - W(d)]^2
##       + γ \sum_n [sum_d G(n,d)q_i - F(n)]^2
##

from dwave.system import LeapHybridSampler
from dwave.system.samplers import DWaveSampler
import dwave_networkx as dnx
import networkx as nx
from dwave.embedding import embed_bqm, embed_qubo, unembed_sampleset
from minorminer import find_embedding
from dimod import BinaryQuadraticModel
from collections import defaultdict
from copy import deepcopy
import pickle

## Setup functions and parameters
## 関数とパラメーターを設定します
### Size parameters
### サイズ パラメーター
numSampling = 1000
for nurses in range(3,5):
    for days in range(5,6):
        #everything below could be a function of `days` and `nurses`
        size = days * nurses

        ### Hard nurse constraint: no nurse on consecutive days
        ### ハード看護師制約：連日出勤は禁止
        a = 7 / 2

        ### Hard shift constraint: enough effort on shift to cover workforce needs
        ### ハード シフト制約：必要なワークフォースを対応できるエフォートの出勤
        lagrange_hard_shift = 1.3
        effort = lambda n : 1.0 # E(n)
        workforce = lambda d : 1.0 # W(d)

        ### Soft nurse constraint: reflect each nurse's preferences
        ### ソフト看護師制約：各々の出勤希望の反映
        lagrange_soft_nurse = 0.3
        preference = lambda n, d : 1.0 # G(n,d)
        duty_days = int(days / nurses) # even distribution

        ### Index function. n = index // days, d = index % days
        ### インデックス関数
        index = lambda n,d: n * days + d

        ## Build Hamiltonian
        ## ハミルトニアンを小袿します
        ### hard nurse constraint: \sum_n,n' \sum_d, d' J_i(n, d)j(n',d')q_iq_j
        ###                       J = a δ_(n,n') δ_(d',d+1)
        J = defaultdict(int)

        for nurse in range(nurses):
            for day in range(days - 1):
                index_d1 = index(nurse, day)
                index_d2 = index(nurse, day + 1)
                J[index_d1, index_d2] = a

        ### Copy to add shift constraints
        ### コピーしてシフトの制約を追加します
        Q = deepcopy(J)

        ### hard shift constraint: λ \sum_d [\sum_n E(n)*q_i - W(d)]^2
        for day in range(days):
            for nurse in range(nurses):
                idx = index(nurse, day)
                Q[idx, idx] += (effort(nurse) - (2 * workforce(day))) * effort(nurse) * lagrange_hard_shift
                for partner in range(nurse +1, nurses):
                    idx2 = index(partner, day)
                    Q[idx, idx2] += 2 * lagrange_hard_shift * effort(nurse) * effort(partner)

        ### soft shift contraint: \sum_n [sum_d G(n,d)q_i - F(n)]^2
        for nurse in range(nurses):
            for day in range(days):
                idx = index(nurse, day)
                Q[idx, idx] += lagrange_soft_nurse * preference(nurse, day) * (preference(nurse, day) - (2 * duty_days))
                for day2 in range(day + 1, days):
                    idx2 = index(nurse, day2)
                    Q[idx, idx2] += 2 * lagrange_soft_nurse * preference(nurse, day) * preference(nurse, day2)


        ## Solve
        ## 解きます

        ### Graph embedding
        chimera_16 = dnx.chimera_graph(16)
        pegasus_16 = dnx.pegasus_graph(4)
        topology = 'pegasus' # 'chimera' or 'pegasus'
        
        fname = "results_%s_N%d_D%d.p" % (topology, nurses, days)
        #previous = pickle.load(open(fname, "rb"))
        #embedding = previous['embedding']
        G = nx.Graph(Q.keys())
        embedding = find_embedding(G.edges, pegasus_16.edges)
        print(G.edges)
        print(embedding)
        print(pegasus_16.edges)
        embeddedQ = embed_qubo(Q, embedding, pegasus_16 )

        ### Energy offset
        ### エネルギー オフセット
        e_offset = lagrange_hard_shift * days * workforce(1) ** 2
        e_offset += lagrange_soft_nurse * nurses * duty_days ** 2

        ### D-Wave sampler
        bqm = BinaryQuadraticModel.from_qubo(embeddedQ, offset=e_offset)
        sbqm = BinaryQuadraticModel.from_qubo(Q, offset=e_offset)

        #hybrid_sampler = LeapHybridSampler()
        sampler = DWaveSampler(solver={'topology__type': topology,'qpu': True})
        results = sampler.sample(bqm, num_reads=numSampling)
        samples = unembed_sampleset(results, embedding, bqm, chain_break_fraction=True)

        ### Save results with pickle for analysis
        fout = "results_%s_N%d_D%d_s%d.p" % (topology, nurses, days, numSampling)
        saveDict = {'results' : results, 'embedding' : embedding, 'bqm': sbqm, 'samples' : samples}
        pickle.dump(saveDict, open(fout, "wb"))

