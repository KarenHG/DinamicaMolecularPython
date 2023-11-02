#subrutina fcc
import numpy as np
rxfcc=np.zeros(5501)
ryfcc=np.zeros(5501)
rzfcc=np.zeros(5501)
rx=np.zeros(5501)
ry=np.zeros(5501)
rz=np.zeros(5501)
vx=np.zeros(5501)
vy=np.zeros(5501)
vz=np.zeros(5501)
rxf=np.zeros(5501)
ryf=np.zeros(5501)
rzf=np.zeros(5501)
n=0

print("crear FCC            [1]")
print("crear foto           [2]")
print("STOP                 [0]")

noption=int(input())

if noption==1:
    print("nc Nat = 4*nc**3")
    nc = int(input())
    nfcc=4*nc**3
    print("número de moléculas",nfcc)
    dens=float(input("densidad "))
    boxx = (nfcc/dens)**(1/3)
    boxy=boxx
    boxz=boxx

    print("boxx",boxx)

    #calcular los lados de la célda unitaria
    cellx=boxx/nc
    celly=boxy/nc
    cellz=boxz/nc
    cell2x=0.5*cellx
    cell2y=0.5*celly
    cell2z=0.5*cellz
    #contruir la celda unitaria
    #subred A
    rxfcc[1]=0.0
    ryfcc[1]=0.0
    rzfcc[1]=0.0
    #subred B
    rxfcc[2]=cell2x
    ryfcc[2]=cell2y
    rzfcc[2]=0.0
    #subred C
    rxfcc[3]=0.0
    ryfcc[3]=cell2y
    rzfcc[3]=cell2z
    #subred D
    rxfcc[4]=cell2x
    ryfcc[4]=0.0
    rzfcc[4]=cell2z
    #construir la red de la celda unitaria
    m=0
    for iz in range(1,nc+1):
        for iy in range(1,nc+1):
            for ix in range(1,nc+1):
                for iref in range(1,5):
                    rxfcc[iref+m]=rxfcc[iref]+cellx*(ix-1)
                    ryfcc[iref+m]=ryfcc[iref]+celly*(iy-1)
                    rzfcc[iref+m]=rzfcc[iref]+cellz*(iz-1)
                m=m+4
                
    #velocidad de las moléculas
    time=0.005
    temp=float(input("Temperatura reducida"))
    print("temperatura reducida ",temp)
    amass=1
    
    for i in range(1,nfcc+1):
        rtemp=np.sqrt(temp/amass)
        vx[i]=rtemp*np.random.normal()
        vy[i]=rtemp*np.random.normal()
        vz[i]=rtemp*np.random.normal()

    #remover el momento neto
    sumx=0
    sumy=0
    sumz=0
    for i in range(1,nfcc+1):
        sumx=sumx+vx[i]*amass
        sumy=sumy+vy[i]*amass
        sumz=sumz+vz[i]*amass
    print("sumx.sumy,sumz",sumx,sumy,sumz)
    print("rtemp",rtemp)
    sumx=sumx/float(nfcc)
    sumy=sumy/float(nfcc)
    sumz=sumz/float(nfcc)
    
    uk=0
    for i in range(1,nfcc+1):
        vx[i]=vx[i]-sumx
        vy[i]=vy[i]-sumy
        vz[i]=vz[i]-sumz
    
    n=nfcc
    for i in range(1,n+1):
        rx[i]=rxfcc[i]
        ry[i]=ryfcc[i]
        rz[i]=rzfcc[i]
    print("crear FCC            [1]")
    print("crear foto           [2]")
    print("STOP                 [0]")

    noption=int(input())

if noption==2:
    file=open("FotoInicial.pdb","w")
    nat=n
    sigmar=3.405
    for i in range(1,n+1):
        rxf[i]=rx[i]*sigmar
        ryf[i]=ry[i]*sigmar
        rzf[i]=rz[i]*sigmar
    for i in range(1,n+1):
        file.write("{:<7} {:<4} {:<7} {:<1} {:<8} {rxf:<7.3f} {ryf:<7.3f} {rzf:<8.3f}\n".format("ATOM",i,"H",1,i,rxf=rxf[i],ryf=ryf[i],rzf=rzf[i]))
    file.close()
    print("crear FCC            [1]")
    print("crear foto           [2]")
    print("STOP                 [0]")

    noption=int(input())

if noption==0:
    name=input("nombre del archivo de salida: ")
    file2=open(name,"w")
    file2.write("{rn:<55} {rbox:<55.50f} {rboxy:<55.50f} \n".format(rn=n,rbox=boxx,rboxy=boxy))
    for i in range(1,n+1):
        file2.write("{rx:<55.50f} {ry:<55.50f} {rz:<55.50f} {vx:<55.50f} {vy:<55.50f} {vz:<55.50f} \n".format(rx=rx[i],ry=ry[i],rz=rz[i],vx=vx[i],vy=vy[i],vz=vz[i]))
    file2.close()
