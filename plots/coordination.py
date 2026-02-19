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

colacn="#D55E00"
colh2o="#2A9D8F"

fig = plt.figure(1, figsize=(4.5, 4.5))
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
    filename = "./results/coordination/%sD%sE%s.cdat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

acnx=pd.DataFrame(widex)
acny=pd.DataFrame(widey)

target="LiPF6inH2O"
widey={}
widex={}
for fi in field:
    filename = "./results/coordination/%sD%sE%s.cdat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

h2ox=pd.DataFrame(widex)
h2oy=pd.DataFrame(widey)


figA=0

if 1:
    kBT=0.592
    colors=sns.color_palette("rocket", 5)
    #acncolors=sns.color_palette("#2A9D8F", n_colors=len(field))
    dark=sns.dark_palette(colors[0], n_colors=1)[0]
    midlight = sns.light_palette(colors[0], n_colors=3)[1]
    acncolors=sns.blend_palette([ midlight, dark], n_colors=len(field))
    acncolors = ramp(colacn, n=len(field), light_to_dark=True)

    
    rmin=3.925
    volume=36*36*36
    density= 511/volume
    coordination=np.zeros((len(acnx.columns)))
    for i, col in enumerate(acny.columns):
        x=acnx[col]
        y=acny[col]
        if(col != "50"):
            plt.plot(x, y, lw=3, color=acncolors[i])
        else:
            plt.plot(x, y, lw=3, color=acncolors[i], label=r"$\mathrm{Li^+ - CH_3CN}$")
        end=78
        coordination[i] = 4*np.pi*density* np.trapezoid(y[:end]*x[:end]*x[:end], x[:end], dx=0.0001)

    #plt.plot(np.arange(len(catx.columns)), coordination)
    
    h2ocolors = ramp(colh2o, n=len(field), light_to_dark=True)

    volume = 28.57**3
    density = 750/volume
    coordinate=np.zeros((len(h2ox.columns)))
    for i, col in enumerate(h2oy.columns):
        x=h2ox[col]
        y=h2oy[col]
        if(col != "50"):
            plt.plot(x,  y, lw=3, color=h2ocolors[i])
        else:
            plt.plot(x, y, lw=3, color=h2ocolors[i], label=r"$\mathrm{Li^+ - H_2O}$")
        end=56
        coordinate[i] = 4*np.pi*density* np.trapezoid(y[:end]*x[:end]*x[:end], x[:end], dx=0.0001)

    rmin=x[57]
 #   plt.plot(x*0+rmin, y, '--')
    plt.yscale("log")    
    plt.xlim(1, 12)
    plt.ylim(0.5e-1, 15)


    plt.xlabel(r'$r \ / \ \mathrm{\AA}$', fontsize=fsize)
    plt.ylabel(r'$g(r)$', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, loc='upper right')
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
#    plt.title(r"$\mathrm{LiPF_6 \ 0.5 \ M}$")
    
#plt.savefig('./figures/rdf.png' , dpi=300, bbox_inches="tight")
plt.show()
