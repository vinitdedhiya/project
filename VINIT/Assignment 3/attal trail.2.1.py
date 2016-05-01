import scipy
# DATA

#1.Henrys constants
Hm=714.286 #(atm/mol)[1]
Hc=29.4 #(atm/mol)[1]

#2.Kla values
Klm=0.004 #(1/s)[2]
Klc=0.004 #(1/s)[2]
Kgw=0.0002 #(1/atm.s)[2]

#3. P_total
P=1.0 #(atm)

#4.Densities
pm=0.656 #(kg/m3)[2]
pc=1.98 #(kg/m3)[2]
pw=1000.0 #(kg/m3)[2]

#5.gas flowrates
Gim=200.0 #(moles)
Gic=200.0 #(moles)
Giw=0.0 #(moles)

#6.liquid flowrates
Lim=500.0 #(moles)
Lic=0.0 #(moles)
Liw=0.0 #(moles)
n=20 #grid points

#7.xinputs
xmi=0
xci=0
xwi=1.0

#8.yinputs
ymi=0.5
yci=0.5
ywi=0

#8.x_guess arrays 
xmguess=n*[xmi]
xcguess=n*[xci]
xwguess=n*[xwi]

#9.y_guess arrays
ymguess=n*[ymi]
ycguess=n*[yci]
ywguess=n*[ywi]

#10.guess array
xyguess=scipy.array(xmguess+xcguess+xwguess+ymguess+ycguess+ywguess)

#11. Diamentions
z=10.0 #(m) 
B=5.0 #(m) 
l=8.0 #(m) 
A=50.0 #(m2)-->z*b

#12.molecular wts
Mm=0.016#kg/mol[2]
Mc=0.044#kg/mol[2]
Mw=0.018#kg/mol[2]

#13. Psat Water at 298K
Pwsat=0.03098

def xm(xc,xw):
    xm=1-xc-xw
   
def xc(xm,xw):
    xc=1-xm-xw
    
def xw(xm,xc):
    xw=1-xm-xc
    
def ym(yc,yw):
    ym=1-yc-yw  
    
def yc(ym,yw):
    yc=1-ym-yw    
    
def yw(ym,yc):
    yw=1-ym-yc    
    
def mix(A,P,Pwsat,Kgw, Mm,Mc,Mw,xm,ym,xc,yc,xw,yw,xmi,ymi,xci,yci,xwi,ywi,Lim,Lic,Liw,Gim,Gic,Giw,pm,pc,pw,Klm,Klc,Hm,Hc):
    n=len(f)
    xm=f[:n/6]
    xc=f[:2*n/6]
    xw=f[:3*n/6]
    ym=f[:4*n/6]
    yc=f[:5*n/6]
    yw=f[5*n/6:]
    dz=z/((n-1))
    ## Forward and backward difference
    errormGLS=Klm*A*((ym[-1]/(Hm))-(xmi*(pm/Mm)))+((ym[-1]-ym[-2])/dz)
    errorcGLS=Klc*A*((yc[-1]/(Hc))-(xci*(pc/Mc)))+((yc[-1]-yc[-2])/dz)
    errorwGLS=((yw[-1]-yw[-2])/dz)-Kgw(yw[-1]-Pwsat)
    errormGUS=Klm*A*((ymi/(Hm))-(xmi*(pm/Mm)))+((ym[1]-ym[0])/dz)
    errorcGUS=Klc*A*((yci/(Hc))-(xci*(pc/Mc)))+((yc[1]-yc[0])/dz)
    errorwGUS=((ywi-yw[0])/dz)-Kgw(ywi-Pwsat) 
    errormLLS=Klm*A*((xm[-1]/(Hm))-(xmi*(pm/Mm)))+((xm[-1]-xm[-2])/dz)
    errorcLLS=Klc*A*((x[-1]/(Hc))-(xci*(pc/Mc)))+((xc[-1]-xc[-2])/dz)
    errorwLLS=((xw[-1]-xw[-2])/dz)-Kgw(xw[-1]-Pwsat)
    errormLUS=Klm*A*((xmi/(Hm))-(xmi*(pm/Mm)))+((xm[1]-xm[0])/dz)
    errorcLUS=Klm*A*((xmi/(Hm))-(xmi*(pm/Mm)))+((xm[1]-xm[0])/dz)
    errorwLUS=((xwi-xw[0])/dz)-Kgw(xwi-Pwsat)
    errormG=scipy.zeros(n/6)
    errormL=scipy.zeros(n/6)
    errorcG=scipy.zeros(n/6)
    errorcL=scipy.zeros(n/6)
    errorwG=scipy.zeros(n/6)
    errorwL=scipy.zeros(n/6)
    errorm[0]=errormGUS;errorm[-1]=errormGLS
    errorc[0]=errorcGUS;errorc[-1]=errorcGLS
    errorw[0]=errorwGUS;errorm[-1]=errorwGLS
    errorm[0]=errormLUS;errorm[-1]=errormLLS
    errorc[0]=errormLUS;errorc[-1]=errorcLLS
    errorw[0]=errorwLUS;errorw[-1]=errorwLLS
    ## Central difference
    errormG[1:-1]=Klm*A*((ym[1:-1]/(Hm))-(xm[1:-1]*(pm/Mm)))+((ym[2:]-ym[1:-1])/dz)
    errorcG[1:-1]=Klc*A*((yc[1:-1]/(Hc))-(xc[1:-1]*(pm/Mm)))+((ym[2:]-ym[1:-1])/dz)
    errorwG[1:-1]=((yw[2:]-yw[1:-1])/dz)-Kgw(yw[1:-1]-Pwsat)
    errormL[1:-1]=Klm*A*((xm[1:-1]/(Hm))-(xm[1:-1]*(pm/Mm)))+((ym[2:]-ym[1:-1])/dz)
    errorcL[1:-1]=Klm*A*((xm[1:-1]/(Hm))-(xm[1:-1]*(pm/Mm)))+((ym[2:]-ym[1:-1])/dz)
    errorwL[1:-1]=((xw[2:]-yw[1:-1])/dz)-Kgw(yw[1:-1]-Pwsat)
    return scipy.concatenate(errormG,errormL,errorcG,errorcL,errorwG,errorwL)

n=len(fguess)
ans=opti.leastsq(mix,fguess,args=(A,P,Pwsat,Kgw, Mm,Mc,Mw,xm,ym,xc,yc,xw,yw,xmi,ymi,xci,yci,xwi,ywi,Lim,Lic,Liw,Gim,Gic,Giw,pm,pc,pw,Klm,Klc,Hm,Hc))
print(ans)
xans=ans[0]
yans=ans[0]
xmans=xm[:n/6]
xmans[0]=xmi
xcans=xc[:2*n/6]
xcans[0]=xci
xwans=xw[:3*n/6]
xwans[0]=xwi
ymans=ym[:4*n/6]
ymans[-1]=ymi
ycans=yc[:5*n/6]
ycans[-1]=yci
ywans=yw[5*n/6:]
ywans[-1]=ywi  
print(xmans)
print(xcans) 
print(xwans) 
print(ymans) 
print(ycans) 
print(ywans)  
    








