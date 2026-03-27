import matplotlib.pyplot as plt

fig, ax = plt.subplots()

start = 1_500_000
plus = 500_000
inflation = 1.05
deposit = 1.1
x,y=[], []
for i in range(20):
    x.append(i+23)
    y.append(start)
    start = (start * (1 + deposit - inflation)) + plus
    # plus = plus * inflation
plt.gca().ticklabel_format(style='plain', axis='y')
ax.plot(x,y)
plt.show()
