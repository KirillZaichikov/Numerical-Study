import math as m

eps = 0.1
alpha = 0.8
beta = 0.2

ksi1 = 0.01
eta1 = 0.1
phi1 = 0.2

print("ksi1, eta1, phi1")
print(ksi1, eta1, phi1)

ksi2 = ksi1 + eps * m.cos(2 * phi1)
eta2 = eta1 + eps * m.sin(2 * phi1)
teta2 = phi1 + 1 + ksi1 + eta1 + eps * m.sin(phi1)

print("ksi2, eta2, teta2")
print(ksi2, eta2, teta2)

##### ОБРАТНОЕ
x = 2.5
# print("da", -teta2 + ksi2 + eta2 + 1)
for i in range(10):
    x_0 = x
    x = x_0 - (x_0+eps*m.sin(x_0)-eps*m.cos(2*x_0)-eps*m.sin(2*x_0)-teta2+ksi2+eta2+1) / \
              (1+eps*m.cos(x_0) + 2*eps*m.sin(2*x_0)-2*eps*m.cos(2*x_0))
    # x = x_0 - (x_0 - eps*m.cos(x_0) - teta2 + ksi2 + eta2 + 1) / (1 + eps*m.sin(x_0))
    print(x)

phi1 = x
ksi1 = ksi2 - eps * m.cos(2*phi1)
eta1 = eta2 - eps * m.sin(2*phi1)

print("ksi1, eta1, phi1")
print(ksi1, eta1, phi1)