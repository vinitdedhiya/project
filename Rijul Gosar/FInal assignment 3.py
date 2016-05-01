ma=1
mb=2
L=10
P=10
U=100
x0=0
y0=400

z1=303

def fof(y,z):
    f=(-P*U*(y-z))/(ma*(4000+0.1*y+0.01*y**2))
    return f
def phi(y,z):
    phi=(-P*U*(y-z))/(mb*(3000+0.2*y+0.05*y**2))
    return phi

tempa=[]
tempb=[]
error=[]
n=[5,6,7,8,9,10,15,30,50]

z0=335
z00=z0
y0=400

    
for j in n:
    y0=400  
    z0=326.5
    h=L/float(j-1)
    print h
    while abs(z1-300)>0.001:           
        for i in range(j):
            tempb.append(z0)            
            k1=h*fof(y0,z0)
            k2=h*fof(y0+k1,z0)
            y0=y0+0.5*(k1+k2)
            k11=h*phi(y0,z0)
            k22=h*phi(y0,z0+k1)
            z0=z0+0.5*(k11+k22)
            tempa.append(y0)        
            tempb.append(z0)
            print y0,z0
    zold=max(tempb)
    z1=z0      
    z0=zold-(z0-300)*0.1
    
    
    
    mcpa=ma*(4000*(400-y0)+0.05*(400**2-y0**2)+(0.01/3)*(400**3-y0**3))
    mcpb=mb*(3000*(z00-z0)+0.1*(z00**2-z0**2)+(0.05/3)*(z00**3-z0**3))
    print mcpa,mcpb
    e=abs(100*(mcpa-mcpb)/mcpa)
    error.append(e)
    
print error
import pylab
pylab.figure("Plot of error vs No. of grids")
pylab.plot(n,error)
pylab.show()
        
    
        
        
            

    

        