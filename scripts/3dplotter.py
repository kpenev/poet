import numpy as np
import matplotlib.pyplot as plt
import sys
class Plotter:
    data=None
    a0_low = a0_high = a0_steps = None
    w0_low = w0_high = w0_steps = None
    def __init__(self, filename):
        f = open(filename)
        line = f.readline()
        split_line = line.split()
        if (split_line[0] != 'a0'):
            print "File is in incorrect format!"
            exit(-1)
        self.a0_low = float(split_line[1])
        self.a0_high = float(split_line[2])
        self.a0_steps = int(split_line[3])
        line = f.readline()
        split_line = line.split()
        if (split_line[0] != 'w0'):
            print "File is in incorrect format!"
            exit(-1)
        self.w0_low = float(split_line[1])
        self.w0_high = float(split_line[2])
        self.w0_steps = int(split_line[3])
        self.data = []
        for i in range(self.a0_steps):
            self.data.append([])
            line = f.readline()
            if line=="": break
            elements = line.split()
            if len(elements) != self.w0_steps:
                print "Line "+str(i)+" "+"has wrong number of elements"
                exit(-1)
            self.data[-1] = [float(x) for x in elements]
        self.data = np.array(self.data)
        self.data[self.data == 100]=np.nan
#        self.data = np.log(self.data)
#        self.data = self.data.clip(0,1000)
#        self.data = 2. - self.data

    def plot(self):
        x = np.linspace(self.w0_low, self.w0_high, self.w0_steps)
        y = np.linspace(self.a0_low, self.a0_high, self.a0_steps)

        Zm = np.ma.masked_where(np.isnan(self.data), self.data)
        plt.pcolormesh(x,y,Zm)
        fontsize = "x-large"
        plt.ylabel("Initial a (AU)", fontsize=fontsize)
        plt.xlabel("Initial $\omega_{conv}$ (rad/day)", fontsize=fontsize)
        plt.title("$log_{10}$ of $L^2$ norm of fractional error",
                  fontsize=fontsize)
        plt.tick_params(labelsize=14)
        cb = plt.colorbar()
        for t in cb.ax.get_yticklabels():
            t.set_fontsize(14)
        plt.show()

if len(sys.argv) != 2:
    print "Usage: " + sys.argv[0]+" filename"
    exit(-1)
p = Plotter(sys.argv[1])
p.plot()
