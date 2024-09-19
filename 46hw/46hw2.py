import numpy as np
import matplotlib.pyplot as plt


t = np.arange(0, np.pi + np.pi/30, np.pi/30)
y = lambda x: np.cos(x)

ty=t*y(t)
s= np.sum(ty)
print('the sum is:',s)


#part b)
import random

# First figure
theta = np.linspace(0, 2 * np.pi, 1000)
R = 1.2
delta_r = 0.1
f = 15
p = 0

x = R * (1 + delta_r * np.sin(f * theta + p)) * np.cos(theta)
y = R * (1 + delta_r * np.sin(f * theta + p)) * np.sin(theta)

plt.figure(figsize=(6, 6))
plt.plot(x, y)
plt.title('Parametric Curve for R = 1.2, δr = 0.1, f = 15, p = 0')
plt.xlabel('x(θ)')
plt.ylabel('y(θ)')
plt.axis('equal')
plt.show()

# Second figure
plt.figure(figsize=(8, 8))
for i in range(1, 11):
    R = i
    delta_r = 0.05
    f = 2 + i
    p = random.uniform(0, 2 * np.pi)
    
    x = R * (1 + delta_r * np.sin(f * theta + p)) * np.cos(theta)
    y = R * (1 + delta_r * np.sin(f * theta + p)) * np.sin(theta)
    
    plt.plot(x, y, label=f'Curve {i}')
    
plt.title('10 Parametric Curves with Varying Parameters')
plt.xlabel('x(θ)')
plt.ylabel('y(θ)')
plt.axis('equal')
plt.legend()
plt.show()
