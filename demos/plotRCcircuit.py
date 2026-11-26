import sys
import numpy as np

prefix = sys.argv[1] if len(sys.argv) > 1 else 'rccircuit'
# columns: time, capacitor voltage U_C, (time state, unused)
data = np.loadtxt('output_test_RCCircuit.txt', usecols=(0, 1))

import matplotlib.pyplot as plt

plt.plot(data[:,0], data[:,1], label='capacitor voltage')
plt.xlabel('time')
plt.ylabel('voltage')
plt.title('RC Circuit Capacitor Voltage')
plt.legend()
plt.grid()
plt.savefig(f'{prefix}_voltage.png', dpi=150)
plt.close()
