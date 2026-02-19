from plotStyle import np, plt, pd, sns
from plotStyle import fsize, lsize, beta, convert, MS, AP, index
from plotStyle import divisionE
import matplotlib.image as mpimg

import matplotlib.colors as mcolors

def ramp(hex, n=5, light_to_dark=True, maxlight=0.75, maxdark=0.):
    base = mcolors.to_rgb(hex)
    def mix(c1, c2, t):
        return tuple((1-t)*a+t*b for a, b in zip(c1, c2))
    white = (1,1,1)
    black = (0,0,0)
    shades=[]
    for i in range(n):
        x=i/(n-1)
        if light_to_dark:
            t_white = (1-x)*maxlight
            t_black = x*maxdark
        else:
            t_white= x*maxlight
            t_black = (1-x)*maxdark
            
        c=mix(base, white, t_white)
        c=mix(c, black, t_black)
        shades.append(mcolors.to_hex(c))
    return shades

fig = plt.figure(1, figsize=(5., 4.5))
plt.tick_params(direction='in', right=True, top=True)
plt.tick_params(labelsize=fsize)
plt.minorticks_on()
plt.tick_params(labelbottom=True, labeltop=False, labelright=False, labelleft=True)
plt.tick_params(direction='in', which='minor', length=5, bottom=True, top=True, right=True)
plt.tick_params(direction='in', which='major', length=10, bottom=True, top=True, right=True)



sigma=2

target="LiPF6inACN"
density="05"
field=["00", "09", "19", "30", "40", "50"]

widey={}
widex={}
for fi in field:
    filename = "../processedData/rdf/%sD%sE%s.dat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

acnx=pd.DataFrame(widex)
acny=pd.DataFrame(widey)

target="LiPF6inH2O"
widey={}
widex={}
for fi in field:
    filename = "../processedData/rdf/%sD%sE%s.dat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

h2ox=pd.DataFrame(widex)
h2oy=pd.DataFrame(widey)

colacn="#D55E00"
colh2o="#2A9D8F"


if 1:
    acncolors = ramp(colacn, n=len(field), light_to_dark=True)
    
    for i, col in enumerate(acny.columns):
        prob = -np.log(acny[col].to_numpy())
        if (col!="50"):
            plt.plot(acnx[col], prob, lw=3, color=acncolors[i])
        else:
            plt.plot(acnx[col], prob, lw=3, color=acncolors[i], label=r"$\mathrm{CH_3CN}$")

    h2ocolors = ramp(colh2o, n=len(field), light_to_dark=True)

    for i, col in enumerate(h2oy.columns):
        prob = -np.log(h2oy[col].to_numpy())
        if (col!="50"):
            plt.plot(h2ox[col], prob, lw=3, color=h2ocolors[i])        
        else:
            plt.plot(h2ox[col], prob, lw=3, color=h2ocolors[i], label=r"$\mathrm{H_2O}$")        


    py = np.linspace(-1, 10, num=100)

    plt.plot(py*0 + 5.3, py, 'k--', lw=2)
    plt.plot(py*0 + 9.1, py, 'k--', lw=2)
        
    plt.xlim(2.4, 11)
    plt.ylim(0.9, 6.75)

    plt.xlabel(r'$R \ / \ \mathrm{\AA}$', fontsize=fsize)
    plt.ylabel(r'$- \ln \rho(R) $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, loc='center')
    
plt.savefig('./figures/fig1c.png' , dpi=300, bbox_inches="tight")
plt.show()
