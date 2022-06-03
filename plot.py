from matplotlib import pyplot as plt
import numpy as np

fig = plt.figure()
fps = 30


tau_max = 5
data = np.loadtxt('./tmp/result_offset_ncc.txt')
data = data.reshape(tau_max, -1, 3)

for t in np.arange(tau_max):
    plt.subplot(5, 2, 2 * t+1)
    plt.plot(data[t][:,0], data[t][:,2], label='tau = {:d}'.format(t+1))
    plt.xlim([-1,5])
    plt.legend()

plt.subplot(2,2,2)
vicon = np.loadtxt('./tmp/vicon_rotation_before.txt')
phone = np.loadtxt('./tmp/phone_rotation_before.txt')

x = np.arange(len(vicon)) / fps
plt.plot(x, vicon, label='vicon')
plt.plot(x, phone, label='phone')
plt.legend()

plt.subplot(2,2,4)
vicon = np.loadtxt('./tmp/vicon_rotation_after.txt')
phone = np.loadtxt('./tmp/phone_rotation_after.txt')

x = np.arange(len(vicon)) / fps
plt.plot(x, vicon, label='vicon')
plt.plot(x, phone, label='phone')
plt.legend()

plt.show()