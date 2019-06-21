import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D

cot = lambda x: 1 / math.tan(x)

def integrate(p, delta, r, theta, phi):
    (omega, kr, ktheta, kphi) = p
    sint = math.sin(theta)
    cost = math.cos(theta)
    cott = cot(theta)
    r_inv = 1 / r
    r_inv_2 = 2 * r_inv
    r_inv_sqr = r_inv * r_inv
    a = 1 - 2 * r_inv
    b = 1 / a

    omega_next  = omega  - delta * (2 * r_inv_sqr * b * omega * kr)
    kr_next     = kr     - delta * (a * r_inv_sqr * omega**2 - b * r_inv_sqr * kr**2 - r * a * ktheta**2 - r * a * sint**2 * kphi**2)
    ktheta_next = ktheta - delta * (r_inv_2 * kr * ktheta - sint * cost * kphi**2)
    kphi_next   = kphi   - delta * (r_inv_2 * kr * kphi + 2 * cott * ktheta * kphi)

    return (omega_next, kr_next, ktheta_next, kphi_next)

count = 1000000
delta = 0.0001

# initial 4-position
r = 2.5
(r, theta, phi) = (r, math.pi / 2, 0)

# initial 4-velocity
omega = 1
kr = 0.075
p = (omega, kr, 0, 1 / r * math.sqrt(((1 - 2 / r)**2 * omega**2 - kr**2) / (1 - 2 / r)))

x = np.ndarray(count)
y = np.ndarray(count)
z = np.ndarray(count)
for i in range(0, count):
    p = integrate(p, delta, r, theta, phi)
    r += delta * p[1]
    theta += delta * p[2]
    phi += delta * p[3]
    
    # renormalize angles
    while theta > math.pi:
        theta -= math.phi
    while phi > 2 * math.pi:
        phi -= 2 * math.pi
    
    # nom
    if r <= 2:
        break

    x[i] = r * math.cos(phi) * math.sin(theta)
    y[i] = r * math.sin(phi) * math.sin(theta)
    z[i] = r * math.cos(theta)
    #light_check = -(1 - 2 / r) * p[0]**2 + p[1]**2 / (1 - 2 / r) + r**2 * p[2]**2 + r**2 * math.sin(theta)**2 * p[3]**2
    if i % 10000 == 0:
        print(i * delta, r, theta, phi)

redshift = omega / p[0]

fig = plt.figure()
ax = fig.gca()
ax.add_artist(plt.Circle((0, 0), 2, color='black'))
ax.plot(x, y, label=f'light path (redshift: {redshift})')
ax.legend()
ax.set_xlim([-13, 13])
ax.set_ylim([-13, 13])
ax.set_aspect('equal')
plt.xlabel('x')
plt.ylabel('y')
plt.show()