import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'verdana'})
rc('text.latex',unicode=True)
rc('text.latex',preamble='sepackage[utf8]{inputenc}')
rc('text.latex',preamble='sepackage[russian]{babel}')


x = np.array([0,	1,	2,	3,	4,	5,	6,	7, 	8,	 9],'float')
y = np.array([29.5,	15,	12,	11,	10,	6.5,	6,	5,	4,	3],'float')

p = y / y.sum()

print(p)

plt.bar(x,p, width=1.0, color="red")
plt.ylabel("состава смеси (доля)")
plt.xlabel("T (C) температура кипения")
plt.show()
