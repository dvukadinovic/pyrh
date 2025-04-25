import pyrh

RH = pyrh.RH(".")

RH.set_keywords()
RH.set_abundances()
RH.set_elements()
RH.dummy()

# print(RH.InputData)
# RH.get_RLK_lines()
# RH.set_InputData()
# print(RH._InputData.Nrays)