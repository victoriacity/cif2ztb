import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt('test_frac.out')
# 'Real' scale factor; gives actual parameters
lj_factor = 20000
# 'Fake' sacle factor; gives scaled parameters
coulomb_factor = 10

dims = [114, 115, 5]

ztb = data[0, :].reshape(dims[0] * dims[1] * dims[2], 3)
ztb_2 = data[1, :].reshape(dims[0] * dims[1] * dims[2], 3)
print(ztb[:10, 2] - np.mean(ztb[:10, 2]))
print(np.min(ztb[:, 2]))

ztb[:, 0] /= lj_factor ** 2
ztb[:, 1] /= lj_factor
ztb[:, 2] /= coulomb_factor
ztb = np.exp(-ztb)
plt.imshow(ztb.reshape(dims[0], dims[1], dims[2], 3, order='C')[:, :, 0, :])
plt.show()
plt.imshow(ztb.reshape(dims[0], dims[1], dims[2], 3, order='C')[:, :, 0, 2])
plt.show()
