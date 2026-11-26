import sys
import numpy as np

prefix = sys.argv[1] if len(sys.argv) > 1 else 'massspring'
data = np.loadtxt('output_test_ode.txt', usecols=(0, 1, 2))

import matplotlib.pyplot as plt

plt.plot(data[:,0], data[:,1], label='position')
plt.plot(data[:,0], data[:,2], label='velocity')
plt.xlabel('time')
plt.ylabel('value')
plt.title('Mass-Spring System Time Evolution')
plt.legend()
plt.grid()
plt.savefig(f'{prefix}_timeevolution.png', dpi=150)
plt.close()

plt.plot(data[:,1], data[:,2], label='phase plot')
plt.xlabel('position')
plt.ylabel('velocity')
plt.title('Mass-Spring System Phase Plot')
plt.legend()
plt.grid()
plt.savefig(f'{prefix}_phase.png', dpi=150)
plt.close()

