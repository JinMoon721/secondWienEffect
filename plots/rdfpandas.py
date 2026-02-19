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
    filename = "./results/rdf/%sD%sE%s.dat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

acnx=pd.DataFrame(widex)
acny=pd.DataFrame(widey)

target="LiPF6inH2O"
widey={}
widex={}
for fi in field:
    filename = "./results/rdf/%sD%sE%s.dat" % (target, density, fi)
    df=pd.read_csv(filename, sep=r"\s+", header=None, names=["r", "hist"])
    widex[fi]= df["r"].to_numpy()
    widey[fi] = df["hist"].to_numpy()

h2ox=pd.DataFrame(widex)
h2oy=pd.DataFrame(widey)

colacn="#D55E00"
colh2o="#2A9D8F"

figA=0

if 1:
    kBT=0.592
    kBT=1
    colors=sns.color_palette("rocket", 5)
    #acncolors=sns.color_palette("#2A9D8F", n_colors=len(field))
    dark=sns.dark_palette(colors[0], n_colors=1)[0]
    midlight = sns.light_palette(colors[0], n_colors=3)[1]
    acncolors=sns.blend_palette([ midlight, dark], n_colors=len(field))
    acncolors = ramp(colacn, n=len(field), light_to_dark=True)
    
    for i, col in enumerate(acny.columns):
        match = np.abs(acnx[col]-5.3).argmin()
        prob = acny[col].to_numpy()
        #prob = acny[col].to_numpy() / acnx[col]**2
        prob -= np.log(prob)
        #prob -= prob[match]
        
        mask = acnx[col] > 5.3
        #prob[mask] -=0.04/kBT
#        plt.plot(acnx[col], 0.1+prob*kBT, lw=3, color=acncolors[i])
        off=0
        if (col!="50"):
            plt.plot(acnx[col], off+prob*kBT, lw=3, color=acncolors[i])
        else:
            plt.plot(acnx[col], off+prob*kBT, lw=3, color=acncolors[i], label=r"$\mathrm{CH_3CN}$")


    h2ocolors=sns.light_palette(colors[2], n_colors=len(field))
    h2ocolors = ramp(colh2o, n=len(field), light_to_dark=True)

    for i, col in enumerate(h2oy.columns):
        match = np.abs(h2ox[col]-5).argmin()
        prob = h2oy[col].to_numpy()
        #prob = h2oy[col].to_numpy()/h2ox[col]**2

        prob -= np.log(prob)
        #prob -= prob[match]
        
        if (col!="50"):
            plt.plot(h2ox[col], prob*kBT, lw=3, color=h2ocolors[i])        
        else:
            plt.plot(h2ox[col], prob*kBT, lw=3, color=h2ocolors[i], label=r"$\mathrm{H_2O}$")        


    py = np.linspace(-0.5, 6, num=100)

    plt.plot(py*0 + 5.3, py/0.592, 'k--', lw=2)
    plt.plot(py*0 + 9.1, py/0.592, 'k--', lw=2)
        
    plt.xlim(2.4, 11)
    plt.ylim(0.6/0.592, 4/0.592)
    #plt.ylim(2, 6)

    plt.xlabel(r'$R \ / \ \mathrm{\AA}$', fontsize=fsize)
    plt.ylabel(r'$- \ln \rho(R) $', fontsize=fsize)
    plt.legend(fontsize=lsize, frameon=False, loc='center')
#    plt.text(-13.5, 3, r"$\textbf{b}$", fontsize=ffsize)
#    plt.title(r"$\mathrm{LiPF_6 \ 0.5 \ M}$")
    
plt.savefig('./figures/rdf.png' , dpi=300, bbox_inches="tight")
plt.show()
