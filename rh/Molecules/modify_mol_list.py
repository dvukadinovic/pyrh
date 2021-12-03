"""
Script to extract data for relevant transitions from large list of lines.

We extract data by min (lmin) and max (lmax) wavelength and the relevant
element isotop in molecule (element_isotop).

In header of file we specify nature of transitions. In UV we care only for
electronic transitions. Following two numbers are Nwavs and qwing used to make
wavelength grid for each molecular line during synthesis (for RH).

Data after 70 character are for lines that are magneticly sensitive. For now,
I have not considered this data since I got an error of "Unapporpriate Hund's
rule" and I am not exactly aware of what is it. And, for now I have not
considered magnetic field in atmosphere, and accordingly I do not need data
for Zeeman splitting of molecular lines.
"""

lines = open("CN/CN_B-X.asc","r").readlines()

out = open("CN/CN_B-X_susi.asc","w")
out_lines = []

lmin, lmax = 300, 410
element_isotop = 12

Nlines = 0
for line in lines:
	line = line.rstrip("\n")[:70]
	lam0 = float(line[:10])
	if lam0>=lmin and lam0<=lmax:
		isotop = int(line[-3:])
		if isotop==element_isotop:
			out_lines.append(line+"\n")
			Nlines += 1

out.writelines(f" {Nlines} MOLECULAR_ELECTRONIC   KURUCZ_NEW\n")
out.writelines("  21    3.25\n")
out.writelines(out_lines)
out.close()