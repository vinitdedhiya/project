ma=1
mb=2
L=10
P=10
U=100
x0=0
y0=400
z0=328.74
z00=z0

def fof(y,z):
    f=(-P*U*(y-z))/(ma*(4000+0.1*y+0.01*y**2))
    return f
def phi(y,z):
    phi=(-P*U*(y-z))/(mb*(3000+0.2*y+0.05*y**2))
    return phi
def cpa(y):
    ca=4000+0.1*y+0.01*y**2
    return ca
def cpb(z):
    cb=3000+0.2*z+0.05*z**2
    return cb
tempa=[]
tempb=[]
error=[]
j_list=[5,6,7,8,9,10,20,30,50,100]
for j in j_list:
    n=j
    print n
    h=L/float(n-1)
    print h    
    
    for i in range(n):
        x=x0+(i*h)
        k1=h*fof(y0,z0) 
        l1=h*phi(y0,z0)
        k2=h*fof(y0+0.5*k1,z0+0.5*l1)
        l2=h*phi(y0+0.5*k1,z0+0.5*l1)
        k3=h*fof(y0+0.5*k2,z0+0.5*l2)
        l3=h*phi(y0+0.5*k2,z0+0.5*l2)
        k4=h*fof(y0+k3,z0+l3)
        l4=h*phi(y0+k3,z0+l3)
        y1=y0+(k1+2*k2+2*k3+k4)/6
        z1=z0+(l1+2*l2+2*l3+l4)/6
        tempa.append(y1)
        tempb.append(z1)
       
        y0=y1
        z0=z1
        #print y1,z1
        tamin=min(tempa)
        tbmax=max(tempb)
        tbmin=min(tempb)
                
        print y1,z1
        print 
    print tamin
    mcpa=ma*(4000*(400-tamin)+0.05*(400**2-tamin**2)+(0.01/3)*(400**3-tamin**3))
    mcpb=mb*(3000*(328.74-tbmin)+0.1*(328.74**2-tbmin**2)+(0.05/3)*(328.74**3-tbmin**3))
    e=100*(mcpa-mcpb)/mcpa
    error.append(e)
    print e
print error
import pylab
pylab.figure("Plot of error vs No. of grids")
pylab.plot(j_list,error)
pylab.show()


    
    
    
    
    
    



