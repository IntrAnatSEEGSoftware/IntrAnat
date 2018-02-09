import math
import pdb
from scipy import ndimage
from numpy import *
from soma import aims



def norme3D(vect):
    return math.sqrt(vect[0]**2+vect[1]**2+vect[2]**2)

def vecteur(a,b):
    return (b[0]-a[0],b[1]-a[1],b[2]-a[2])

def findMoyMax(npCT,sizex,sizey,sizez,volCT,entry,d,s1n,coteRegionX,coteRegionY,coteRegionZ,CT):
    #initialization of variables
    maxi=0
    inc=0
    if CT==True:
        compX1=s1n[0]/sizex-coteRegionX
        compX2=s1n[0]/sizex+coteRegionX
        compY1=s1n[1]/sizey-coteRegionY
        compY2=s1n[1]/sizey+coteRegionY
        compZ1=s1n[2]/sizez-coteRegionZ
        compZ2=s1n[2]/sizez+coteRegionZ
        if compZ1 <=0:
            compZ1=0
        if compY1 <=0:
            compY1=0
        if compX1 <=0:
            compX1=0    
        if compX2>= len(npCT[0][0])-1:
            compX2= len(npCT[0][0])-1  
        if compY2>= len(npCT[0][0])-1:
            compY2= len(npCT[0][0])-1  
        if compZ2>= len(npCT[0][0])-1:
            compZ2= len(npCT[0][0])-1  
    else:
        compX1=s1n[0]-coteRegionX
        compX2=s1n[0]+coteRegionX
        compY1=s1n[1]-coteRegionY
        compY2=s1n[1]+coteRegionY
        compZ1=s1n[2]-coteRegionZ
        compZ2=s1n[2]+coteRegionZ
        if compZ1 <=0:
            compZ1=0
        if compY1 <=0:
            compY1=0
        if compX1 <=0:
            compX1=0    
        if compX2>= len(npCT[0][0])-1:
            compX2= len(npCT[0][0])-1  
        if compY2>= len(npCT[0][0])-1:
            compY2= len(npCT[0][0])-1  
        if compZ2>= len(npCT[0][0])-1:
            compZ2= len(npCT[0][0])-1     
            
        print "compX1",compX1,"compX2",compX2,"compY1",compY1,"compY2",compY2,"compZ1",compZ1,"compZ2",compZ2    
            
    #we will do the computation until we find a value higher than 1500
    while maxi==0 and inc<1000:
        #print "compX1",compX1,"compX2",compX2,"compY1",compY1,"compY2",compY2,"compZ1",compZ1,"compZ2",compZ2
        #computation in a CT
        if CT==True:
            #we take out any values less than 2500
            newnpCT=npCT[0,round(compZ1):round(compZ2),round(compY1):round(compY2),round(compX1):round(compX2)]
            newnpCT=newnpCT.clip(1500)
            newnpCT[newnpCT==1500]=0
            #We make an opening
            CTopened=ndimage.grey_opening(newnpCT, size=(int(round(1.4/sizez)),int(round(1.4/sizey)),int(round(1.4/sizex)))) #en mm ca soit isotropic
            newnpCT=npCT[0,round(compZ1):round(compZ2),round(compY1):round(compY2),round(compX1):round(compX2)]
            newnpCT=newnpCT.clip(1500)
            newnpCT[newnpCT==1500]=0
        else:
            #we take out any values above 50
            newnpCT=npCT[0,round(compZ1):round(compZ2),round(compY1):round(compY2),round(compX1):round(compX2)]
            newnpCT=newnpCT.clip(0,50)
            newnpCT[newnpCT==50]=0
            #We make an opening
            CTopened=ndimage.grey_opening(newnpCT, size=(int(round(1.4/sizez)),int(round(1.4/sizey)),int(round(1.4/sizex)))) #en mm ca soit isotropic
            newnpCT=npCT[0,round(compZ1):round(compZ2),round(compY1):round(compY2),round(compX1):round(compX2)]
            newnpCT=newnpCT.clip(0,100)
            newnpCT[newnpCT==100]=0
        try:
          maxi=CTopened.max()
        except:
          pdb.set_trace()
        #if the opening was too strong and didn't keep any values but the original array had value>1500 we keep the newnpCT (original one)  
        if maxi==0 and newnpCT.max()!=0:
            CTopened=newnpCT
            maxi=newnpCT.max()
        #if there is nothing near we open the ROI    
        if (maxi<2500 and CT==True) or (CT==False and maxi>50 ) or maxi==0:
            compX1-=sizex
            if compX1<= 0:
                compX1= 0
                
            compX2+=sizex
            if compX2>= len(npCT[0][0][0])-1:
                compX2= len(npCT[0][0][0])-1
                
            compY1-=sizey
            if compY1<= 0:
                compY1= 0

            compY2+=sizey    
            if compY2>= len(npCT[0][0])-1:
                compY2= len(npCT[0][0])-1              
            
            compZ1-=sizez
            if compZ1<= 0:
                compZ1= 0    
            
            compZ2+=sizez 
            if compZ2>= len(npCT[0])-1:
                compZ2= len(npCT[0])-1
            
            maxi=0
                 
        inc+=1       
        
    #calculation of the center of mass
    (lbl,numfeatures)=ndimage.label(CTopened)
    centre=ndimage.measurements.center_of_mass(CTopened,lbl,range(1,numfeatures+1))
    #if there is one center
    if numfeatures==1:
        moy = ndimage.measurements.center_of_mass(CTopened)
    #if there are more than 5 centers of mass we put moy to none    
    elif numfeatures>5:
        print "trop de centres de masse"
        moy=None
    #if there are between 3 and 5 centers of mass we take the 2 bigger ones    
    elif numfeatures>2 and numfeatures<6:
        p=1
        to={}
        suppr=[]
        #we determine the size of each centers
        while p<=lbl.max():
            val=where(lbl==p)
            tailleVal={p:len(val[0])}
            to.update(tailleVal)
            p+=1
        #sort the sizes and only keep the 2 biggest ones    
        tailles=[x for x in to.values()]
        tailles.sort()
        del tailles[-2:]
        
        #then we know wich centers to suppress...
        for el in tailles:
            for cle, valeur in tailleVal.items():
                if valeur==el:
                    suppr.append(cle)
        #...and suppress them
        for el in suppr:
                lbl[lbl==el]=0
        #we then compute the centers of mass, since we do not know on wich axis they are, we first try with [0,1], if one of the centers is (nan,nan,nan) we try with [1,2]
        centre=ndimage.measurements.center_of_mass(CTopened,lbl,[0,1])
        if isnan(centre[0][0])==True:
            centre=ndimage.measurements.center_of_mass(CTopened,lbl,[1,2])
            
        #we calculate the coordinates in the CT native space in order to calculate the distance with the first approximation
        centre0=(centre[0][2]+compX1*sizex,centre[0][1]+compY1*sizey,centre[0][0]+compZ1*sizez)
        centre1=(centre[1][2]+compX1*sizex,centre[1][1]+compY1*sizey,centre[1][0]+compZ1*sizez)
        entCentre0=vecteur(entry,centre0)
        dist0=norme3D(entCentre0)
        entCentre1=vecteur(entry,centre1)
        dist1=norme3D(entCentre1)
        #we only take the one that is nearest to the point we first approximated
        if dist0>dist1:
            moy=centre[1]
        else:
            moy=centre[0]
    #here is the case where nothing goes wrong        
    else:
        centre0=(centre[0][2]+compX1*sizex,centre[0][1]+compY1*sizey,centre[0][0]+compZ1*sizez)
        centre1=(centre[1][2]+compX1*sizex,centre[1][1]+compY1*sizey,centre[1][0]+compZ1*sizez)
        entCentre0=vecteur(entry,centre0)
        dist0=norme3D(entCentre0)
        entCentre1=vecteur(entry,centre1)
        dist1=norme3D(entCentre1)
        if dist0>dist1:
            moy=centre[1]
        else:
            moy=centre[0]
    #if we have a center we will also calculate its distance with the previous plot       
    if moy is not None:        
        moyverif=((moy[2]+compX1)*sizex,(moy[1]+compY1)*sizey,(moy[0]+compZ1)*sizez)   
        verif=norme3D(vecteur(entry,moyverif))   
        emVect=vecteur(entry,moyverif)
        #if the plot found is too far we put moy to none
        if d!=0 and verif>d*1.3:
            moy=None
            print "too far"
        #this is the case where we approximate the target, we want it at most 2mm far from the theoretical target
        elif d==0 and verif>2:
            moy=None
            print "verif:" ,verif
        #if the plot is too close     
        elif d!=0 and verif<d/1.3:
            print "trop pres"
            #pdb.set_trace()
            #coef2=d*1.1/verif
            #moy=((coef2*(moyverif[2]-entry[2])+entry[2])/sizex-compX1,(coef2*(moyverif[1]-entry[1])+entry[1])/sizey-compY1,(coef2*(moyverif[0]-entry[0])+entry[0])/sizez-compZ1)
            moy=None
            
    if moy is not None:
            moy0=moy[0]
            moy1=moy[1]
            moy2=moy[2]
            moy0+=compZ1
            moy1+=compY1
            moy2+=compX1
            moy=(moy0,moy1,moy2)
            
    return moy

def locateContact(npCT,coteRegionX,coteRegionY,coteRegionZ,sizex,sizey,sizez,volCT,entry,d,point,s1n,serpentin,CT):
        #initialization of variables
        error=10   
        j=0
        moyt=(0,0,0)
        ite=0
        #print "entry: ", entry
        #print "point : ",point
        #we do the center of mass of the ROI while the error between two approximations is < 0.04
        while error>0.04 and ite<10000:  
            #initialization of variables
            npCTtempo=None 
            aa=0
            
            #We launch the approximation a first time in order to know if moy is none
            if j==0:
                moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,entry,d,s1n,coteRegionX,coteRegionY,coteRegionZ,CT)
            else:
                moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,entry,d,replace,coteRegionX,coteRegionY,coteRegionZ,CT)
              
            #print "moy",moy
            #if moy is none we want to get out of the while
            if moy is None:
                error=0
            else:    
                error=abs(norme3D(moyt)-norme3D(moy))
                replace=(moy[2]*sizex,moy[1]*sizey,moy[0]*sizez)
                moyt=moy
                ite+=1
            j+=1    
        #we lauch the approximation of the center of mass once         
        if moy is None:
            if s1n==point:
                moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,point,0,s1n,coteRegionX,coteRegionY,coteRegionZ,CT)
            else:    
                moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,point,0.2,s1n,coteRegionX,coteRegionY,coteRegionZ,CT)
            #if the computations fail at this point, we will give it the first approximation
            if moy is None:
                appPointret=s1n
                        
            #transformation of the center of mass found to the CT natif space    
            else:            
                appPointret=(moy[2]*sizex,moy[1]*sizey,moy[0]*sizez)  
        else:           
            appPointret=(moy[2]*sizex,moy[1]*sizey,moy[0]*sizez)  
        #print "appPointret: ",appPointret    
        return appPointret

    
def locateContacts(target,entry,npCT,volCT,nbContacts,sizex,sizey,sizez,do,transfo_pre_to_postopInv,brainMask,sizeT1,dicPoints,serpentin,transfo_pre_to_postop,CT):
    #entry et target dans le repere CT natif
    #variables' initialization
    targetH=target
    entryH=entry
    contacts={}
    i=0
    theta=0
    signex={}
    signey={}
    signez={}
    angles={}
    theta=0
    #counts the number of times moy is returned none
    compteurMoy=0
    
    #Approximation for each contacts, we start at the theoretical target, wich will also be approximated
    #it has to be noted that the entry is after the first iteration the current approximated plot, and target the previous one
    while i<nbContacts:
        print "contact numero:", i
        #first and second approximation are different from others
        #the first entry in the do dictionnary, wich rassembles the inter-contact's distances, is the length between the target and the entry, we then have to take the distance
        #between the target and the next contact: do[1]
        if i==0:
            coteRegionX=abs(do[i+1]/(1.7*sizex))
            coteRegionY=abs(do[i+1]/(1.7*sizey))
            coteRegionZ=abs(do[i+1]/(1.7*sizez))
        #after the first approximation we just take the inter contact distance with the point we want to approximate and the previous one    
        else:
            coteRegionX=abs(do[i]/(1.7*sizex))
            coteRegionY=abs(do[i]/(1.7*sizey)) 
            coteRegionZ=abs(do[i]/(1.7*sizez))
              
        #at first, we take the theoretical target as the first estimation    
        if i==0:
            s1n=target
        #otherwise we calculate the next point with distance d from the previous plot    
        else:    
            if i==1:
                a=(do[1]/do[0])
                #our first approximation in the CT natif space
                s1n=(a*(entry[0]-target[0])+target[0],a*(entry[1]-target[1])+target[1],a*(entry[2]-target[2])+target[2])
            else:
                a=(do[i]/do[i-1])
                #our first approximation in the CT natif space
                s1n=(a*(entry[0]-target[0])+entry[0],a*(entry[1]-target[1])+entry[1],a*(entry[2]-target[2])+entry[2])      
        if i>1:
            aa=1.0005
            b=0.9994
            enS1=vecteur(entry,s1n)
            enS1Norm=norme3D(enS1)
            while enS1Norm<(do[i]-0.4) or enS1Norm>(do[i]+0.4):
                if enS1Norm<do[i]-0.4:
                    s1n=((a*(entry[0]-target[0])+entry[0])*aa,(a*(entry[1]-target[1])+entry[1])*aa,(a*(entry[2]-target[2])+entry[2])*aa)
                    aa+=0.0001 
                    enS1=vecteur(entry,s1n)
                    enS1Norm=norme3D(enS1) 
                else:
                    s1n=((a*(entry[0]-target[0])+entry[0])*b,(a*(entry[1]-target[1])+entry[1])*b,(a*(entry[2]-target[2])+entry[2])*b)
                    b-=0.0001
                    enS1=vecteur(entry,s1n)
                    enS1Norm=norme3D(enS1) 
                    if b<=0:
                        s1n=(2*entry[0]-target[0],2*entry[1]-target[1],2*entry[2]-target[2])
                        enS1Norm=do[i]
                   
        #print "dicPoints: ",dicPoints           
        #print "s1n :", s1n
        #first approximation with d=0, distance between target and its approximation should be small
        #returns a value in the CT natif space without voxel size correction
        if i==0:
            appPoint=locateContact(npCT,coteRegionX,coteRegionY,coteRegionZ,sizex,sizey,sizez,volCT,target,0,dicPoints[i+1],s1n,serpentin,CT)
        #second opproximation, target is the previous point, it will be entry after this iteration    
        elif i==1:
            appPoint=locateContact(npCT,coteRegionX,coteRegionY,coteRegionZ,sizex,sizey,sizez,volCT,target,do[i],dicPoints[i+1],s1n,serpentin,CT)
        #entry is the previous point    
        else:
            appPoint=locateContact(npCT,coteRegionX,coteRegionY,coteRegionZ,sizex,sizey,sizez,volCT,entry,do[i],dicPoints[i+1],s1n,serpentin,CT)
        #at the second iteration the target becomes the entry (wich is from iteration 2 the previous point)    
        if i>1:
            prepre=target
            target=entry
            #print "target 312 : ", target
        if i>0:    
            entry=appPoint 
            #print "entry 315: ",entry
            vect11=vecteur(target,entry)
        #Calculation of the angle between the two vectors joining 3 consecutive points, last one being the current approximated point.    
        if i>1:
            try:
                theta=math.acos(vdot(vect11,vect12)/(norme3D(vect11)*norme3D(vect12))) 
                #print theta
            #if it can't be done it is often because the vectors are identical, so we instanciate theta to 0    
            except:
                theta=0    
                

        #We store the angles in order to be able to reduce deviations            
        if i>1:   
            signex.update({i:((appPoint[0]-target[0])/abs(appPoint[0]-target[0]))})
            signey.update({i:((appPoint[1]-target[1])/abs(appPoint[1]-target[1]))})
            signez.update({i:((appPoint[2]-target[2])/abs(appPoint[2]-target[2]))})
            angle={i:theta}
            angles.update(angle)
                
        #target becomes the approximated one    
        if i==0:
            target=appPoint
            #print "target 338: ", target
        if i>0:
            vect12=vect11      
            entry=appPoint 
            #print "entry 342: ",entry
            
        #print entry    
        #Transformation of the found plot to the T1 natif referential  
        appPointtemp=(appPoint[0],appPoint[1],appPoint[2])
        appPointtemp=list(appPointtemp)
        appPointtemp.append(1)
        appPointtemp=array(appPointtemp)
        appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
        point=list(appPointT1nat)
        del point[-1]
        appPointT1nat=tuple(point)
        contact={i:appPointT1nat}
        contacts.update(contact)
        
        if serpentin==True:
            if i>2:
                if (signex[i-1]-signex[i]!=0 or signey[i-1]-signey[i]!=0 or signey[i-1]-signey[i]!=0) and angles[i]>0.05 and angles[i-1]>0.05:
                    prepreNorm=norme3D(vecteur(prepre,entry))
                    coef=do[i-1]/prepreNorm                   
                    target=((entry[0]-prepre[0])*coef+prepre[0],(entry[1]-prepre[1])*coef+prepre[1],(entry[2]-prepre[2])*coef+prepre[2])
                    #print "target 363 : ", target
                    appPointtemp=(target[0],target[1],target[2])
                    appPointtemp=list(appPointtemp)
                    appPointtemp.append(1)
                    appPointtemp=array(appPointtemp)
                    appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
                    point=list(appPointT1nat)
                    del point[-1]
                    appPointT1nat=tuple(point)
                    contacts[i-1]=appPointT1nat
                    print "modif 410"
                    
        else:
            if appPoint==s1n:
                #increases the number of times moy is none
                compteurMoy+=1
                if compteurMoy>1:
                    print "trop de moy=None"   
                    moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,dicPoints[i+1],0.2,dicPoints[i+1],coteRegionX,coteRegionY,coteRegionZ,CT)    
                    if moy is not None:  
                        appPoint=(moy[2]*sizex,moy[1]*sizey,moy[0]*sizez) 
                        print "appPoint: ", appPoint
                        #print "entry 383: ", entry
                    else:
                        appPoint=dicPoints[i+1]
                        #print "entry 386 :",entry
                    compteurMoy=0
                    #Transformation of the found plot to the T1 natif referential  
                    appPointtemp=(appPoint[0],appPoint[1],appPoint[2])
                    appPointtemp=list(appPointtemp)
                    appPointtemp.append(1)
                    appPointtemp=array(appPointtemp)
                    appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
                    point=list(appPointT1nat)
                    del point[-1]
                    appPointT1nat=tuple(point)
                    contacts[i]=appPointT1nat
                    print "modif 437"
                    try:
                        coef1=(do[i-1]+do[i-2])/norme3D(vecteur(contacts[i-3],contacts[i]))
                        contacts[i-1]=(coef1*(appPointT1nat[0]-contacts[i-3][0])+contacts[i-3][0],coef1*(appPointT1nat[1]-contacts[i-3][1])+contacts[i-3][1],coef1*(appPointT1nat[2]-contacts[i-3][2])+contacts[i-3][2])
                        coef2=do[i-2]/norme3D(vecteur(contacts[i-3],contacts[i]))
                        contacts[i-2]=(coef2*(appPointT1nat[0]-contacts[i-3][0])+contacts[i-3][0],coef2*(appPointT1nat[1]-contacts[i-3][1])+contacts[i-3][1],coef2*(appPointT1nat[2]-contacts[i-3][2])+contacts[i-3][2])    
                    except:
                        try:
                            coef2=do[i-1]/norme3D(vecteur(contacts[i-2],contacts[i]))
                            contacts[i-1]=(coef2*(appPointT1nat[0]-contacts[i-2][0])+contacts[i-2][0],coef2*(appPointT1nat[1]-contacts[i-2][1])+contacts[i-2][1],coef2*(appPointT1nat[2]-contacts[i-2][2])+contacts[i-2][2])   
                        except:
                            print "modif 450"
                    #print "a"    
                    entry=appPoint
                    #print "entry 408: ", entry
                    #Transformation of the found plot to the T1 natif referential  
                    appPointtemp=(contacts[i-1][0],contacts[i-1][1],contacts[i-1][2])
                    appPointtemp=list(appPointtemp)
                    appPointtemp.append(1)
                    appPointtemp=array(appPointtemp)
                    appPointT1nat=transfo_pre_to_postop.dot(appPointtemp.T)
                    point=list(appPointT1nat)
                    del point[-1]
                    target=tuple(point)
                    appPointT1nat=contacts[i]
                    #print "target 418: ", target
                    
        #we make shure that the new point isn't at an angle>10 degrees
        if i>=nbContacts-2:
            v0=vecteur(contacts[i-2],contacts[i-1])
            v1=vecteur(contacts[i-1],contacts[i])
            try:
                theta=math.acos(vdot(v0,v1)/(norme3D(v0)*norme3D(v1))) 
            except:
                theta=0
            if theta> 0.174533:
                #we lauch the approximation of the center of mass once          
                moy=findMoyMax(npCT,sizex,sizey,sizez,volCT,entry,do[i],dicPoints[i+1],coteRegionX/1.1,coteRegionY/1.1,coteRegionZ/1.1,CT)
                #if the computations fail at this point, we will give it the first approximation
                if moy is None:
                    appPoint=dicPoints[i+1]
                else:
                    appPoint=(moy[2]*sizex,moy[1]*sizey,moy[0]*sizez)     
                #Transformation of the approximated contact to the T1 natif referential
                appPointtemp=list(appPoint)
                appPointtemp.append(1)
                appPointtemp=array(appPointtemp)
                appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
                point=list(appPointT1nat)
                del point[-1]
                contacts[i]=tuple(point)  
                print "modif 490"
                v0=vecteur(contacts[i-2],contacts[i-1])
                v1=vecteur(contacts[i-1],contacts[i])
                try:
                    theta=math.acos(vdot(v0,v1)/(norme3D(v0)*norme3D(v1))) 
                except:
                    theta=0
                if theta> 0.174533:
                    appPoint=dicPoints[i+1]
                    c#Transformation of the approximated contact to the T1 natif referential
                    appPointtemp=list(appPoint)
                    appPointtemp.append(1)
                    appPointtemp=array(appPointtemp)
                    appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
                    point=list(appPointT1nat)
                    del point[-1]
                    contacts[i]=tuple(point)     
                    print "modif 507"
                #print "theta sup 0.17"
            
        #print "target 445: ", target
        #print "entry :", entry
        #Suppression of the located contact so it doesn't interfear with the approximation of the next one
        #Calculation of the size of the npCT matrix we are going to remove
        if i==0:
            a=round(appPoint[2]/sizez-do[i+1]/(2.2*sizez))
            b=round(appPoint[2]/sizez+do[i+1]/(2.2*sizez))
            c=round(appPoint[1]/sizey-do[i+1]/(2.2*sizey))
            d=round(appPoint[1]/sizey+do[i+1]/(2.2*sizey))
            e=round(appPoint[0]/sizex-do[i+1]/(2.2*sizex))
            f=round(appPoint[0]/sizex+do[i+1]/(2.2*sizex))
        else:    
            a=round(appPoint[2]/sizez-do[i]/(2.2*sizez))
            b=round(appPoint[2]/sizez+do[i]/(2.2*sizez))
            c=round(appPoint[1]/sizey-do[i]/(2.2*sizey))
            d=round(appPoint[1]/sizey+do[i]/(2.2*sizey))
            e=round(appPoint[0]/sizex-do[i]/(2.2*sizex))
            f=round(appPoint[0]/sizex+do[i]/(2.2*sizex))
        #Removal of the contact        
        try:
            npCT[0,a:b,c:d,e:f]=zeros(npCT[0,a:b,c:d,e:f].shape)
        except:
            pdb.set_trace()
        print "contacts : ", contacts    
        
        #approximation of the contacts with only the first approximation :newS1=(a*(entry[0]-target[0])+entry[0],a*(entry[1]-target[1])+entry[1],a*(entry[2]-target[2])+entry[2]) if the contact is near the bone
        #Only starts running after the nb of contacts/2-th- contact
        if i>nbContacts/2:
            #if no brainMask is found, this computation doesn't take place
            if brainMask is None:
                pass
            else:   
                #We see if we still are in the brain: the value will be !=0
                if brainMask[0,appPointT1nat[2],appPointT1nat[1],appPointT1nat[0]]!=0:
                    pass
                else:
                    #if we weren't in brain, we are going to look around if we find some
                    if brainMask[0,appPointT1nat[2]-(do[i]/sizeT1[2]):appPointT1nat[2]+(do[i]/sizeT1[2]),appPointT1nat[1]-(do[i]/sizeT1[1]):appPointT1nat[1]+(do[i]/sizeT1[1]),appPointT1nat[0]-(do[i]/sizeT1[0]):appPointT1nat[0]+(do[i]/sizeT1[0])].max()!=0:
                        pass
                    else: 
                        #no brain is found, the contact will then be approximated in the continuation of the 2 contacts that come before
                        i+=1
                        while i<nbContacts:
                            #print "b"
                            entar=vecteur(target,entry)
                            entarNorm=norme3D(entar)
                            a=(do[i]/entarNorm)
                            newS1=(a*(entry[0]-target[0])+entry[0],a*(entry[1]-target[1])+entry[1],a*(entry[2]-target[2])+entry[2]) 
                            
                            #Transformation of the approximated contact to the T1 natif referential
                            appPointtemp=list(newS1)
                            appPointtemp.append(1)
                            appPointtemp=array(appPointtemp)
                            appPointT1nat=transfo_pre_to_postopInv.dot(appPointtemp.T)
                            point=list(appPointT1nat)
                            del point[-1]
                            newContact=tuple(point)            
                            contact={i:newContact}
                            contacts.update(contact)
                            target=entry
                            #print "target 504: ", target
                            entry=newS1
                            #print "entry 506: ",entry
                            i+=1
        i+=1   
    return contacts