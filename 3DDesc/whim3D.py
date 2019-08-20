import scipy
import scipy.linalg
from .AtomProperty import GetRelativeAtomicProperty
from .Atom3DProperty import GetAtomCoordinateMatrix

Version=1.0
#############################################################################




def XPreCenter(X):
    """
    #################################################################
    Center the data matrix X
    #################################################################
    """
    Xdim=scipy.size(X,axis=0)
    Xmean=scipy.mean(X,axis=0)
    Xmean=scipy.matrix(Xmean)
    Xp=X-scipy.ones([Xdim,1])*Xmean
   
    return Xp


def GetPropertyMatrix(AtomLabel,proname='m'):
    """
    #################################################################
    #################################################################
    """
    res=[]
    for i in AtomLabel:
        res.append(GetRelativeAtomicProperty(i,proname))
    
    return scipy.matrix(scipy.diag(res))
        


def GetSVDEig(CoordinateMatrix,AtomLabel,proname='u'):
    """
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    
    S=XPreCenter(CoordinateMatrix)

    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
    
    return s




def GetWHIM1(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    --->L1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0],3)


def GetWHIM2(CoordinateMatrix,AtomLabel,proname='u'):
    """
    #################################################################
    WHIM descriptors
    
    --->L2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1],3)


def GetWHIM3(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->L3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[2],3)


def GetWHIM4(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Tu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    T=round(sum(s),3)
    
    return T

def GetWHIM5(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Au
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    
    return round(A,3)


def GetWHIM6(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Vu
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    A=s[0]*s[1]+s[0]*s[2]+s[1]*s[2]
    T=sum(s)
    V=A+T+s[0]*s[1]*s[2]
    
    return round(V,3)

def GetWHIM7(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P1u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[0]/(s[0]+s[1]+s[2]),3)


def GetWHIM8(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->P2u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    return round(s[1]/(s[0]+s[1]+s[2]),3)




def GetWHIM9(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Ku
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)
    
    res=0.0
    for i in s:
        res=res+abs(i/sum(s)-1/3.0)
        
    Ku=3.0/4*res
    
    return round(Ku,3) 


def GetWHIM10(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E1u
    #################################################################
    """
    
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[0],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,0]).T,4))
    
    return round(float(res.real),3)
    
    
def GetWHIM11(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E2u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    res=scipy.power(s[1],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,1]).T,4))
    
    return round(float(res.real),3)
    

def GetWHIM12(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->E3u
    #################################################################
    """
    nAtom,kc= CoordinateMatrix.shape   
    
    if proname=='u':
        weight=scipy.matrix(scipy.eye(nAtom))
    else:      
        weight=GetPropertyMatrix(AtomLabel,proname)
    S=XPreCenter(CoordinateMatrix)
    u,s,v=scipy.linalg.svd(S.T*weight*S/sum(scipy.diag(weight)))
     
    try:
        res=scipy.power(s[2],2)*nAtom/sum(scipy.power(S*scipy.matrix(u[:,2]).T,4))
        return round(float(res.real),3)
    except: return "NA"
    
def GetWHIM13(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors
    
    --->Du
    #################################################################
    """
    c1=GetWHIM10(CoordinateMatrix,AtomLabel,proname)
    c2=GetWHIM11(CoordinateMatrix,AtomLabel,proname)
    c3=GetWHIM12(CoordinateMatrix,AtomLabel,proname)
    Du=c1+c2+c3
    
    return round(float(Du),3)




def GetWHIM14(CoordinateMatrix,AtomLabel,proname='u'):
    
    """
    #################################################################
    WHIM descriptors

    --->P3u
    #################################################################
    """
    s=GetSVDEig(CoordinateMatrix,AtomLabel,proname)

    return round(s[2]/(s[0]+s[1]+s[2]),3)



def GetWhim3D(lcoordinates):

    dout = {}
    CoordinateMatrix, AtomLabel = GetAtomCoordinateMatrix(lcoordinates)

    dout['L1u'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='u')
    dout['L2u'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='u')
    dout['L3u'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='u')
    dout['Tu'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='u')
    dout['Au'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='u')
    dout['Vu'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='u')
    dout['P1u'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='u')
    dout['P2u'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='u')
    dout['Ku'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='u')
    dout['E1u'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='u')
    dout['E2u'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='u')
    dout['E3u'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='u')
    dout['Du'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='u')
    dout['L1m'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='m')
    dout['L2m'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='m')
    dout['L3m'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='m')
    dout['Tm'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='m')
    dout['Am'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='m')
    dout['Vm'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='m')
    dout['P1m'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='m')
    dout['P2m'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='m')
    dout['Km'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='m')
    dout['E1m'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='m')
    dout['E2m'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='m')
    dout['E3m'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='m')
    dout['Dm'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='m')
    dout['L1e'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='En')
    dout['L2e'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='En')
    dout['L3e'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='En')
    dout['Te'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='En')
    dout['Ae'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='En')
    dout['Ve'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='En')
    dout['P1e'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='En')
    dout['P2e'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='En')
    dout['Ke'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='En')
    dout['E1e'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='En')
    dout['E2e'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='En')
    dout['E3e'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='En')
    dout['De'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='En')
    dout['L1v'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='V')
    try:
        dout['L2v'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='V')
    except:
        dout['L2v'] = "NA"
    dout['L3v'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='V')
    dout['Tv'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='V')
    dout['Av'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='V')
    dout['Vv'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='V')
    dout['P1v'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='V')
    dout['P2v'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='V')
    dout['Kv'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='V')
    dout['E1v'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='V')
    dout['E2v'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='V')
    dout['E3v'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='V')
    dout['Dv'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='V')
    dout['L1p'] = GetWHIM1(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['L2p'] = GetWHIM2(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['L3p'] = GetWHIM3(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['Tp'] = GetWHIM4(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['Ap'] = GetWHIM5(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['Vp'] = GetWHIM6(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['P1p'] = GetWHIM7(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['P2p'] = GetWHIM8(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['Kp'] = GetWHIM9(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['E1p'] = GetWHIM10(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['E2p'] = GetWHIM11(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['E3p'] = GetWHIM12(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['Dp'] = GetWHIM13(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['P3p'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='alapha')
    dout['P3u'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='u')
    dout['P3m'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='m')
    dout['P3e'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='En')
    dout['P3v'] = GetWHIM14(CoordinateMatrix, AtomLabel, proname='V')

    return dout