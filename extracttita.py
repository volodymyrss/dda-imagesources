import pyfits
from numpy import *
import pywcs
import sys
from fitimage import FitMosaicSources
import os


mosaic=pyfits.open(sys.argv[1])

wcs=pywcs.WCS(mosaic[4].header)

class Image:

    all_ebins=None
    inclusive=False

    flux=None
    variance=None

    def __init__(self,e1,e2):
        self.e1=e1
        self.e2=e2
        self.all_ebins=[]

    def __repr__(self):
        return "[Image %.5lg %.5lg]"%(self.e1,self.e2)

    def add(self,e1,e2,flux,variance,exposure,header):
        print repr(self),"requested to add",e1,e2

        if not ((not self.inclusive and (self.e1<=e1<=self.e2 and self.e1<=e2<=self.e2)) or \
                (self.inclusive and (self.e1<=e1<=self.e2 or self.e1<=e2<=self.e2))):
            return False

        if self.flux is None:
            print repr(self),"init"
            self.flux=flux
            self.variance=variance
            self.exposure=exposure
            self.header=header
        else:
            print repr(self),"adding"
            self.flux+=flux
            self.variance+=variance

        self.significance=self.flux/self.variance**0.5
        self.all_ebins.append([e1,e2])
    
    def get_mebins(self):
        return min(zip(*self.all_ebins)[0]),max(zip(*self.all_ebins)[1])

    def write(self):
        hdus=[pyfits.PrimaryHDU(),pyfits.ImageHDU(),pyfits.ImageHDU(self.flux),pyfits.ImageHDU(self.variance),pyfits.ImageHDU(self.significance),pyfits.ImageHDU(self.exposure)]
        

        cts=hdus[2].data

        print "arf in this bin",self.arf_bin
        flux=cts/self.arf_bin
        flux_hdu=pyfits.ImageHDU(flux)
        hdus.append(flux_hdu)
        
        me1,me2=m.get_mebins()
        for h in hdus[2:]:
            h.header=self.header
            h.header['E_MIN']=me1
            h.header['E_MAX']=me2

        pyfits.HDUList(hdus).writeto("%s_%.5lg_%.5lg.fits"%(name.replace(" ","_"),me1,me2),clobber=True)



if __name__=="__main__":
    arf_e1=pyfits.open("/Integral/data/resources/arfs/arf_62_1528_11_pre.fits")[1].data['ENERG_LO']
    arf_e2=pyfits.open("/Integral/data/resources/arfs/arf_62_1528_11_pre.fits")[1].data['ENERG_HI']
    #arf=pyfits.open("/Integral/data/resources/arfs/arf_62_1528_11_pre.fits")[1].data['SPECRESP']
    arf=pyfits.open("/Integral/data/resources/arfs/arf_mc_oldenergies.fits")[1].data['SPECRESP']
    from scipy.interpolate import interp1d

    for l in open(os.path.dirname(os.path.realpath(__file__))+"/remnants.txt"):
        ra,dec,name=l.strip().split(",")
        ra=float(ra)
        dec=float(dec)

        x,y=wcs.wcs_sky2pix([ra],[dec],0)
        x=x[0]
        y=y[0]

        sx,sy=mosaic[4].data.shape 

        if not ( 0<x<sx and 0<y<sy ): 
            continue

        print name,ra,dec,x,y
        print "in image"
        
        re1=[mosaic[2+4*ie].header['E_MIN'] for ie in range((len(mosaic)-2)/4)]
        re2=[mosaic[2+4*ie].header['E_MAX'] for ie in range((len(mosaic)-2)/4)]

        merged=[]
       # merged=[[(48,69.5),None],[(69.5,75.5),None],[(75.5,110),None]]
        merged+=[Image(a,b) for a,b in zip(re1,re2)]
        #merged=[[(48,65.5),None],[(65.5,69.5),None],[(69.5,73.5),None],[(73.5,75.5),None],[(75.5,79.5),None],[(79.5,85.5),None],[(82,110),None]]
        merged+=[Image(65.5,85.5),Image(20,50)]
 #       merged+=[[(24,48),None],[(48,65.5),None],[(65.5,85.5),None],[(85.5,109.5),None]]


        for ie in range((len(mosaic)-2)/4):
            e1=mosaic[2+4*ie].header['E_MIN']
            e2=mosaic[2+4*ie].header['E_MAX']
            print "image",e1,e2

            for m in merged:
                m.add(e1,e2,copy(mosaic[2+4*ie].data),copy(mosaic[3+4*ie].data),copy(mosaic[5+4*ie].data),header=mosaic[2+4*ie].header.copy())
            

        fluxes=[]
        for m in merged:
            print "fluxes for",repr(m),m.get_mebins(),m.all_ebins
            me1,me2=m.get_mebins()
            if m.flux is None: continue

            fl,vr,ex=m.flux,m.variance,m.exposure
            ex=mosaic[-1].data

            cx,cy=meshgrid(arange(fl.shape[0]),arange(fl.shape[1]))

            cp=copy(fl)
            cp[(x-cx)**2+(y-cy)**2>3**2]=-1000

            print x,cx.max(),cp.max(),cp.argmax(),cp.shape
            peak=unravel_index(cp.argmax(),fl.shape)
            print "peak at",peak,x,y,round(x),round(y)
            print "peak vs x,y,",fl[peak[0],peak[1]],fl[round(x),round(y)]
            print "around",fl[262,262],fl[261,261]
            flux=fl[peak]
            error=vr[peak]**0.5
            flux2=fl[round(x),round(y)]
            error2=vr[x,y]**0.5


     #       flux3,error3=extract_flux(fl,vr,x,y)
            fit=FitMosaicSources()
            fit.free_constant=False
            fit.e1=me1
            fit.e2=me2
            fit.wcs=pywcs.WCS(m.header)
            #fit.cat=[[name,(peak[0],peak[1])]]
            fit.cat=[[name,(ra,dec)]]
            #fit.usexy=True
            r=fit.fit_image(fl,vr)
            flux3=r[1][0]
            error3=r[2][0]

            constant=0

            fit.free_constant=True
            r=fit.fit_image(fl,vr)
            flux4=r[1][0]
            error4=r[2][0]

            print me1,me2,flux,error
    


            mec=(me1+me2)/2.
            dme=(me2-me1) #/2.

            m.arf_bin=interp1d(arf_e1,arf)(mec)
            m.write()

            fluxes.append([mec,dme,flux/dme,error/dme,flux/error,ex[x,y],flux2/dme,error2/dme,flux3/dme,error3/dme,flux4/dme,error4/dme,constant,peak[0],peak[1],peak[0]-y,peak[1]-x,m.arf_bin])

        savetxt('%s_fluxes.txt'%name.replace(" ","_"),fluxes,fmt="%.5lg")

        reg="""
# Region file format: DS9 version 4.1
# Filename: significance63.5_65.5.fits
global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
fk5
circle(%.5lg,%.5lg,1080") # text = {%s}
"""
        open("%s.reg"%name,"w").write(reg%(ra,dec,name))



