import pyrh
import pickle

t = pyrh.Test()
print(t.get_again())
send_me = pickle.dumps(t)
del t
new = pickle.loads(send_me)
print(new.get_again())

# RH = pyrh.RH(".")

# a = pickle.dumps(RH.input)

# RH.set_keywords()
# RH.set_abundances()
# RH.set_elements()
# RH.dummy()