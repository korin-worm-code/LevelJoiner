good = []
bad = False
for r in range(0,9):
	value = stdev[r]
	if value >= 25:
		bad = True
		continue
	else:
		if bad == True:
			good = []
			bad = False
		good += [value]
print good