import pyrh
from schwimmbad import MPIPool

def parallelise_me(RH):
    # RH.set_keywords()
    # RH.set_abundances()
    # RH.set_elements()
    # RH.dummy()
    print(RH.get(2))

with MPIPool() as pool:
    # RH = pyrh.RH(".")
    RH = pyrh.Test()

    pool.map(parallelise_me, [RH]*10)