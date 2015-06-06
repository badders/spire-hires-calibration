from pylab import *
from scipy.constants import h, k, c

rcParams['font.family'] = 'serif'
rcParams['font.sans-serif'] = 'CMU Sans Serif'
rcParams['font.serif'] = 'CMU Serif'
rcParams['font.monospace'] = 'CMU Typewriter'

Ts = [2.73, 10, 100, 300, 5500]

def bb(l):
    return (2 * h * c**2 / l**5) * (1 /(e**(h*c / (l * k * T)) - 1))

def bb_f(f, T):
    return (2 * h * f**3 / c**2) * (1 /(e**(h*f / (k * T)) - 1))

l = logspace(7, 16, num=1000)

for T in Ts:
    loglog(l, bb_f(l, T), label='{}K'.format(T))

ylabel('Intensity / $W.m^{-2}Hz^{-1}sr^{-1}$')
xlabel('Frequency / $Hz$')
ylim(1e-25, 1e-7)
legend(loc='best')
grid()

savefig('./figures/bb10k.pdf')

figure()
ylabel('Intensity / $W.m^{-2}Hz^{-1}sr^{-1}$')
xlabel('Frequency / $Hz$')

wls = array([70, 100, 160, 250, 350, 500]) * 1e-6
fs = c / wls

l = logspace(10, 13, num=1000)

loglog(l, bb_f(l, 10), 'b-')
loglog(fs, bb_f(fs, 10), 'ro', markersize=6)

ylim(1e-25, 1e-15)
grid()

savefig('./figures/bb10kHERSCHEL.pdf')

show()
