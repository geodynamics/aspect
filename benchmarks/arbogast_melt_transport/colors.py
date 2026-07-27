#dark green
#blue
#pink
#red
#orange
#light green
c = ['#0A5F02', '#1c567a', '#814292', '#d7383b', '#fdae61', '#c0f8b8']
m = ['o', 'D', 's','>','<','^','v']
def marker(idx):
    return m[idx%len(m)]

def color(idx):
    return c[idx%len(c)]


if __name__ == "__main__":

    import os, sys, numpy as np, matplotlib.pyplot as plt

    x=np.array([0, 100])
    y1=np.array([1,1])
    plt.plot(x,y1,'+-', color=color(0), label='color0', linewidth=3.0)
    plt.plot(x,y1+1,'+-', color=color(1), label='color1', linewidth=3.0)
    plt.plot(x,y1+2,'+-', color=color(2), label='color2', linewidth=3.0)
    plt.plot(x,y1+3,'+-', color=color(3), label='color3', linewidth=3.0)
    plt.plot(x,y1+4,'+-', color=color(4), label='color4', linewidth=3.0)
    plt.plot(x,y1+5,'+-', color=color(5), label='color5', linewidth=3.0)

    plt.ylim([0,10])
    plt.legend()
    plt.show()
