import numpy as np
import matplotlib.pyplot as plt
import numba 

#parametros
maxnat=5000
maxlist=9000000
maxnesp=3
maxnpair=6
maxpairs=6
deltar=0.05
maxg=90000
maxk=900000
gr=np.zeros(10000)
r0=np.zeros(10000)

#leemos param.mix
natesp=np.zeros(3,dtype=int)
sigma=np.zeros((3,3))
eps=np.zeros((3,3))
amass=np.zeros(maxnat)
amassi=np.zeros(3)
iesp=np.zeros(maxnat)
file = open("param.mix","r")
nesp=int(file.readline())
sum_nat=0.0
for i in range(1,nesp+1):
    itype=file.readline()
    partes = file.readline().split()
    natesp[i]=int(partes[0])
    amassi[i]=float(partes[1])
    sigma[i,i]=float(partes[2])
    eps[i,i]=float(partes[3])
    sum_nat=sum_nat+natesp[i]
file.close() 
nat=int(sum_nat)
npp=0
for i in range(1,nesp+1):
    for j in range(i,nesp+1):
        npp = npp+1
        sigma[i,j] = (sigma[i,i]+sigma[j,j])/2
        eps[i,j] = np.sqrt(eps[i,i]*eps[j,j])
        eps[j,i] = eps[i,j]
        sigma[j,i] = sigma[i,j]
        
print("\n i   j   sigma   epsilon")
for i in range(1,nesp+1):
    for j in range(i,nesp+1):
        print("\n "+str(i)+"   "+str(j)+"    "+str(round(sigma[i,j],3))+"    "+str(round(eps[i,j],3)))
  
n=0
for i in range(1,nesp+1):
    for j in range(1,int(natesp[i]+1)):
        n=n+1
        iesp[n]=i
        amass[n]=amassi[i]      
        
#calcular la fracción molar inicial

xmol=np.zeros(nesp+1)
for i in range(1,nesp+1):
    xmol[i]=natesp[i]/sum_nat
    print("\nnmol y xmol de la especie ",i,natesp[i],xmol[i])
    
#leer el archivo md.dat(parámetros de entrada)

file2 = open("md.dat","r")
cnfile_original=file2.readline()
cnfile=cnfile_original.replace("\n","")
nstep=int(file2.readline())
nprint=int(file2.readline())
ngr=int(file2.readline())
nfilm=int(file2.readline())
time=float(file2.readline())
temp=float(file2.readline())
rcut=float(file2.readline())

# Leemos la configuración inicial

file = open(cnfile,"r")
a=file.readline().split()
n=int(a[0])
box=boxx=boxy=boxz=float(a[1])

rx=np.zeros(n+1)
ry=np.zeros(n+1)
rz=np.zeros(n+1)
vx=np.zeros(n+1)
vy=np.zeros(n+1)
vz=np.zeros(n+1)

for i in range(1,n+1):
    partes = file.readline().split()
    rx[i]=float(partes[0])
    ry[i]=float(partes[1])
    rz[i]=float(partes[2])
    vx[i]=float(partes[3])
    vy[i]=float(partes[4])
    vz[i]=float(partes[5])
file.close()

boxix=boxiy=boxiz=1/box

#variables iniciales en cero

nga=np.zeros(maxlist)      
pnorm=np.zeros(maxg)
ptang=np.zeros(maxg)

temptotal=0.0
ekitotal=0.0
epitotal=0.0
pressit=0.0
etotalt=0.0

temptotal2=0.0
ekitotal2=0.0
epitotal2=0.0
pressit2=0.0
etotal2=0.0

#correcciones de largo alcance

def lrc(nat,dens,box):
    #for i in range(1,nesp+1):
    #   sigmaSuma=sigmaSuma+natesp[i]*sigma[i,i]
    #   sigmaPromedio=sigmaSuma/nat
    sigmaChecar=1.0 #lo metemos dentro de cclos para trabajar cada sigma?
    sr3 = (sigmaChecar/rcut)**3
    sr9 = sr3**3
    Ulrc12 = 8*np.pi*dens*float(nat)*sr9/9
    Ulrc6 = -8*np.pi*dens*float(nat)*sr3/3
    Ulrc = Ulrc12 + Ulrc6
    plrc12 = 4*Ulrc12
    plrc6 = 2*Ulrc6
    plrc = plrc12 + plrc6
    
    Ulrc=Ulrc/nat
    vol = box**3
    plrc= plrc/vol
    
    print("Correcciones de largo alcance")
    print("Ulrc: ",round(Ulrc,4))
    print("plrc: ",round(plrc,4))
    return(Ulrc,plrc)

#imprimir los valores de entrada

print("\nPotencial de Lennard-Jones\n")
print("Dinámica molecular NVT contante\n")

vol=boxx*boxy*boxz
rho=dens=nat/vol

print("Número de paso de tiempo             ",nstep)
print("Frecuencia de impresión de salida    ",nprint)
print("Número de moléculas                  ",nat)
print("Densidad reducida                    ",rho)
print("Tiempo de pasos reducido             ",time)
print("Temperatura deseada                  ",temp)
print("Distancia de radion de corte         ",rcut)
print("Longitud de la caja                  ",boxx,boxy,boxz)
Ulrc,plrc=lrc(nat,dens,box)
print("\n   ISTEP   TEMPI  PRESSI  EKI    EPI    ET  ")

fx=np.zeros(nat+1)
fy=np.zeros(nat+1)
fz=np.zeros(nat+1)
@numba.jit(nopython=True)
def fuerza(rx,ry,rz):
    fx=np.zeros(nat+1)
    fy=np.zeros(nat+1)
    fz=np.zeros(nat+1)
    deltaz = 0.05
    rcutsq = rcut**2
    vir=0.0
    epi=0.0
    virxx=0.0
    viryy=0.0
    virzz=0.0
    
    for i in range(1,nat):
        for j in range(i+1,nat+1):
            
            iei=int(iesp[i])
            iej=int(iesp[j])
            
            dx = rx[i]-rx[j]
            dy = ry[i]-ry[j]
            dz = rz[i]-rz[j]
            
            dx = dx-np.round(dx*boxix)*boxx
            dy = dy-np.round(dy*boxiy)*boxy
            dz = dz-np.round(dz*boxiz)*boxz
            
            rijsq= dx**2 + dy**2 + dz**2
            rij = np.sqrt(rijsq)
            
            if(rij<rcut):
                s6=(sigma[iei,iej]/rij)**6
                uij=4*eps[iei,iej]*s6*(s6-1)
                epi=epi+uij
                duij=-24*eps[iei,iej]*s6*(2*s6-1)/rij
                
                fijx = -duij*dx/rij
                fijy = -duij*dy/rij
                fijz = -duij*dz/rij
                
                virxx = virxx + dx*fijx
                viryy = viryy + dy*fijy
                virzz = virzz + dz*fijz
                
                fx[j] = fx[j] - fijx
                fy[j] = fy[j] - fijy
                fz[j] = fz[j] - fijz
                
                fx[i] = fx[i] + fijx
                fy[i] = fy[i] + fijy
                fz[i] = fz[i] + fijz
          
    return(fx,fy,fz,epi,virxx,viryy,virzz)

vxt=np.zeros(nat+1)
vyt=np.zeros(nat+1)
vzt=np.zeros(nat+1)            

def mover(istepp):
    global rx
    global ry
    global rz
    global vx
    global vy
    global vz
    
    deltaz=0.05
    timesq=time**2
               
    wxxkk=0.0
    wyykk=0.0
    wzzkk=0.0
    eki=0.0
    zside=-0.5*boxz
    factf=np.zeros(nat+1)
    vxx=np.zeros(nat+1)
    vyy=np.zeros(nat+1)
    vzz=np.zeros(nat+1)
    
    
    #valores de la posición al tiempo t
    
    vxt=vx
    vyt=vy
    vzt=vz
        
    factf[1:nat+1]=time/amass[1:nat+1]
    vxx[1:nat+1] = vx[1:nat+1] + 0.5*time*(fx[1:nat+1]/amass[1:nat+1])
    vyy[1:nat+1] = vy[1:nat+1] + 0.5*time*(fy[1:nat+1]/amass[1:nat+1])
    vzz[1:nat+1] = vz[1:nat+1] + 0.5*time*(fz[1:nat+1]/amass[1:nat+1])    
        
        #predecir velocidades, posiciones y calcular la energía cinética
        
    vx = vx + fx*factf
    vy = vy + fy*factf
    vz = vz + fz*factf
        
    rx = rx + vx*time
    ry = ry + vy*time
    rz = rz + vz*time
        
    rx = rx - np.round(rx*boxix)*boxx
    ry = ry - np.round(ry*boxiy)*boxy
    rz = rz - np.round(rz*boxiz)*boxz
        
    eki=0.0
    
    for i in range(1,nat+1):
    
        vvx = 0.5*(vxt[i]+vx[i])
        vvy = 0.5*(vyt[i]+vy[i])
        vvz = 0.5*(vzt[i]+vz[i])
        eki = eki + amass[i]*(vvx**2+vvy**2+vvz**2)
        
        wxxkk = wxxkk + amass[i]*vvx**2
        wyykk = wyykk + amass[i]*vvy**2
        wzzkk = wzzkk + amass[i]*vvz**2
        
        vtang = 0.5*amass[i]*(vvx**2+vvy**2)
        vnorm = amass[i]*vvz**2
        
        nbinz = int(abs(rz[i]-zside)/deltaz+1)
        pnorm[nbinz] = pnorm[nbinz] +vnorm
        ptang[nbinz] = ptang[nbinz] +vtang
        
    wvxx = wxxkk
    wvyy = wyykk
    wvzz = wzzkk
    
    eki = 0.5*eki
    dof = 3*nat -3
    tempi = 2*eki/dof
    vol = boxx*boxy*boxz
    
    #reescalar las velocidades para matener la temperatura constante
    
    #remover el momento neto   
    sumx=0.0
    sumy=0.0
    sumz=0.0
    sumtot=0.0
    
    
    sumx = np.sum(vx*amass[0:nat+1])
    sumy = np.sum(vy*amass[0:nat+1])
    sumz = np.sum(vz*amass[0:nat+1])
    sumtot = np.sum(amass)
        
    sumx=sumx/sumtot
    sumy=sumy/sumtot
    sumz=sumz/sumtot
    
    
    vx=vx-sumx
    vy=vy-sumy
    vz=vz-sumz
    
    ratio = np.sqrt(temp/tempi)
    
    
    vx = vx*ratio
    vy = vy*ratio
    vz = vz*ratio
      
    return(eki,tempi,wvxx,wvyy,wvzz)

def instant(istepp,tempii,ekii,epii):
    global temptotal
    global ekitotal
    global epitotal
    global pressit
    global etotalt
    global temptotal2
    global ekitotal2
    global epitotal2
    global pressit2
    global etotal2 
    vol=boxx*boxy*boxz
    
    pxxt = virxx + wvxx
    pyyt = viryy + wvyy
    pzzt = virzz + wvzz
    
    pressi = (pxxt +pyyt + pzzt)/3/vol
    
    epii = epii/nat    
    ekii = ekii/nat
    etotal = ekii + epii
    
    if(istepp%nprint==0):
        print("  ",istepp,"   ",round(tempii,4),round(pressi,4),round(ekii,4),round(epii,4),round(etotal,4))
      
    temptotal = temptotal  + tempii
    ekitotal  = ekitotal   + ekii
    epitotal  = epitotal   + epii
    pressit   = pressit  + pressi  
    etotalt   = etotalt    + etotal
                                                                        
    temptotal2= temptotal2 + tempii**2
    ekitotal2 = ekitotal2  + ekii**2
    epitotal2 = epitotal2  + epii**2
    pressit2  = pressit2  + pressit**2 
    etotal2   = etotal2    + etotal**2
    
    return(pressi,ekii,epii,etotal)
    
def promedios():
    global dens
    vol=boxx*boxy*boxz
    dens=nat/vol
    
    rmean_t = temptotal/nstep
    rmean_eki = ekitotal/nstep
    rmean_epi = epitotal/nstep + Ulrc
    rmean_pressi = pressit/nstep + plrc
    rmean_etotal = etotalt/nstep
    rmean_epi_sin_corre = epitotal/nstep
    rmean_pressi_sin_corre = pressit/nstep
    
    fluct_t = np.sqrt(abs(temptotal2/nstep - rmean_t**2))
    fluct_eki = np.sqrt(abs(ekitotal2/nstep - rmean_eki**2))
    fluct_epi = np.sqrt(abs(epitotal2/nstep - rmean_epi**2))
    fluct_pressi = np.sqrt(abs(pressit2/nstep - rmean_pressi**2))
    fluct_etot = np.sqrt(abs(etotal2/nstep - rmean_etotal**2))
    
    print("\n               PROMEDIOS               ")
    print("<T>  =",round(rmean_t,5),"+-",round(fluct_t,6))
    print("<K>  =",round(rmean_eki,5),"+-",round(fluct_eki,6))
    print("<EP> =",round(rmean_epi,5),"+-",round(fluct_epi,6))
    print("<ET> =",round(rmean_etotal,5),"+-",round(fluct_etot,6))
    print("<Pres>=",round(rmean_pressi,5),"+-",round(fluct_pressi,6))
    print(round(rmean_epi_sin_corre,5),round(rmean_pressi_sin_corre,5))
    print("\n Densidad",round(dens,4))
    
#########################################################
for istep in range(1,nstep+1):
    fx,fy,fz,epi,virxx,viryy,virzz=fuerza(rx,ry,rz)
    eki,tempi,wvxx,wvyy,wvzz=mover(istep) 
    pressi,eki,epi,etotal=instant(istep,tempi,eki,epi)
promedios()
