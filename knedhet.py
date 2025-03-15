import matplotlib.pyplot as plt
import numpy as np
import time
import distinctipy

mu = np.linspace(0, 1, 1000)
nu = np.linspace(0.5, 1, 1000)
A = 0.95
x_ = 0.0
calc_time = 1000
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(nu[0], nu[-1])
ax.set_ylim(mu[0], mu[-1])

def map(x, mu, nu):
    return -mu+A*abs(x)**nu

x, y, categories= [], [], []
example_line1 = '01'*int(calc_time / 2)
example_line2 = '10'*int(calc_time / 2)
for i, par1 in enumerate(nu):
    print(i)
    for j, par2 in enumerate(mu):
        times = 1
        x_=0
        counter = 0
        while True:
            tmp_line=''
            for k in range(calc_time):
                x_old = x_
                for count in range(times):
                    x_ = map(x_, par2, par1)
                if x_ - x_old > 0:
                    tmp_line += '1'
                else:
                    tmp_line += '0'
            if (tmp_line != example_line1 and tmp_line != example_line2) or times > 1024:
                x.append(par1)
                y.append(par2)
                categories.append(counter)
                break
            else:
                times *= 2
                counter += 1

with open('heteroclinics\\'+'A'+f'{A}'[:3]+'nu'+f'{nu[0]}'[:3]+f'{nu[-1]}'[:3]+'.txt', 'w') as dat:
    dat.write(str(nu[0])) # nu
    dat.write('\n')
    dat.write(str(nu[-1])) # nu
    dat.write('\n')
    dat.write(str(mu[0])) # nu
    dat.write('\n')
    dat.write(str(mu[-1])) # nu
    dat.write('\n')
    dat.write(str(len(mu))) # nu
    dat.write('\n')
    for cat in categories:
        dat.write(str(cat))
        dat.write('\n')

print(set(categories))
colormap = np.array([(0.29703107260951633, 0.2957743936235912, 0.9272974416894829), (0.3183314280884725, 0.9323252362960627, 0.3610405053601042),
                     (0.9659505112472502, 0.35957237433430556, 0.2935091699454408), (0.9423834820124043, 0.9002023970023137, 0.28710211598563307), 
                     (0.29261556185970256, 0.9452968001374702, 0.9752137765023777),  (0.3306538384641366, 0.37968777604908516, 0.29385792924882637), 
                     (0.9655097820641895, 0.3171010391554942, 0.9968303145994776), (0.620791207342849, 0.6359034791441474, 0.7081490128276622), 
                     (0.6420537420076042, 0.2883177473871127, 0.5933514940330531), (0.713775692504449, 0.9891023976464999, 0.6651576931812739),
                     (0.6928939747103388, 0.6037612457478428, 0.2919065043660857), (0.9478100549301239, 0.6791273156989974, 0.9548547714343467), 
                     (0.3061949467333815, 0.639765467265412, 0.9846404555922371), (0.29054476098665255, 0.6065480673108111, 0.5389376170300988), 
                     (0.9997096455531786, 0.6980911124599046, 0.554805744120247), (0.6916324488545669, 0.4560221933703202, 0.9618466866261876), 
                     (0.3518161059842052, 0.8591562321968021, 0.6775081421518184), (0.6399608892687929, 0.8845942056981821, 0.9837447207260359), 
                     (0.3009723020629114, 0.29370718586430605, 0.5923007494767217), (0.6392564675766446, 0.991849494721768, 0.32537849789849593)])

ax.scatter(x, y, c=colormap[categories], s=1)
plt.show()