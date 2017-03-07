import astropy.io.fits as pyfits
from numpy import *
import astropy.wcs as pywcs
import sys

class FitMosaicSources:

    #cat=[('Tycho SNR',(6.3583,64.1528))] 
    cat=[('Cas A',(350.85,58.810))] 

    radius=5
    free_constant=False

    usexy=False

    #cat=[('Cen B',(206.704330,-60.391499)),('4U 1344-60',(206.883331,-60.610001))] 

    def model(self):
        if not hasattr(self,'wcs'):
            raise

       # print "data shape",self.data_shape

        model=zeros(self.data_shape) # x,y?..
        x,y=meshgrid(arange(model.shape[0]),arange(model.shape[1]))
        sigma=self.sigma

        center=None

        self.centers=[]

      #  print self.cat,self.fluxes

        for (sn,(ra,dec)),flux in zip(self.cat,self.fluxes):
       #     print sn,ra,dec,flux

            if self.usexy:
                sx,sy=[ra],[dec]
            else:
                sx,sy=self.wcs.wcs_world2pix([ra],[dec],0)
        #    print sx[0],sy[0]

            self.centers.append([sx[0],sy[0]])

            if center is None:
                center={'x':sx[0],'y':sy[0],'n':1}
            else:
                center['n']+=1
                center['x']=(center['x']*(center['n']-1)+sx[0])/center['n']
                center['y']=(center['y']*(center['n']-1)+sy[0])/center['n']

            model+=exp(-((x-sx)**2+(y-sy)**2)/sigma**2/2)*flux

        if self.free_constant:
            model+=self.constant

        c=center['x'],center['y'],self.radius
        #print "center",c
        #mask=(abs(x-c[0])<c[2]) & (abs(y-c[1])<c[2])
        mask=((x-c[0])**2+(y-c[1])**2)<c[2]**2

        return model,mask


    def fit_image(self,rate,variance):
        self.fluxes=ones((len(self.cat),))
        self.data_shape=rate.shape

        self.sigma=1.1578793
        self.free_sigma=False

        import nlopt
            

        def myfunc(x, grad):
            #print ":::::::",x

            if self.free_sigma:
                self.sigma=x[0]
                x=x[1:]

            if self.free_constant:
                self.constant=x[0]
                x=x[1:]

            self.fluxes=x[:]

            model,mask=self.model()
            r=(rate-model)/variance**0.5
            self.residuals=r
            self.residuals[~mask]=0
            ndof=r[mask].flatten().shape[0]
            residual=(r[mask]**2).sum()/ndof
            #print "------------->",x,model.sum(),residual,ndof,residual/ndof
            return residual



        #self.free_constant=True

        x0=[]
        xmin=[]
        xmax=[]

        if self.free_sigma:
            x0.append(1)
            xmin.append(0.1)
            xmax.append(2)
        
        if self.free_constant:
            x0.append(0)
            xmin.append(-10)
            xmax.append(10)

        x0+=[0]*len(self.fluxes)
        xmin+=[-100]*len(self.fluxes)
        xmax+=[1000]*len(self.fluxes)

        opt = nlopt.opt(nlopt.LN_BOBYQA, len(x0))
        #opt = nlopt.opt(nlopt.LN_COBYLA, len(x0))
        opt.set_lower_bounds(xmin)
        opt.set_upper_bounds(xmax)
        opt.set_min_objective(myfunc)
        opt.set_xtol_rel(1e-4)
        x = opt.optimize(x0)

        if self.free_sigma:
            self.sigma=x[0]
            x=x[1:]

        if self.free_constant:
            self.constant=x[0]
            x=x[1:]

        self.fluxes=x[:]

        print self.centers
        #print variance
        for c in self.centers:
            print c,c[0],c[1],variance[int(c[0]),int(c[1])]**0.5 

        self.flux_errors=[variance[int(c[0]),int(c[1])]**0.5 for c in self.centers]

        #print c,variance[c[0],c[1]].shape

        minf = opt.last_optimum_value()
        print "optimum at ", x
        print "minimum value = ", minf
        print "result code = ", opt.last_optimize_result()

        model=self.model()[0]

        if False:
            pyfits.HDUList([
                        pyfits.PrimaryHDU(rate-model,header=e.header),
                        pyfits.ImageHDU(model,header=e.header),
                        pyfits.ImageHDU(rate,header=e.header),
                        pyfits.ImageHDU(variance,header=e.header)
                        ]).writeto("residuals_%.5lg_%.5lg.fits"%(self.e1,self.e2),clobber=True)

        return [self.e1,self.e2],self.fluxes,self.flux_errors

    def main(self):
        f=pyfits.open(self.mosaic_fn)
    
        intensity_i={}
        variance_i={}
        for ie,e in enumerate(f):
            if 'IMATYPE' in e.header and e.header['IMATYPE']=="INTENSITY":
                self.e1,self.e2=[e.header[a] for a in ['E_MIN','E_MAX']]
                intensity_i[(self.e1,self.e2)]=ie
            
            if 'IMATYPE' in e.header and e.header['IMATYPE']=="VARIANCE":
                self.e1,self.e2=[e.header[a] for a in ['E_MIN','E_MAX']]
                variance_i[(self.e1,self.e2)]=ie
                
        spectra=[] 

        for self.e1,self.e2 in intensity_i.keys():
            print self.e1,self.e2
            e=f[intensity_i[(self.e1,self.e2)]]
            rate=e.data
            variance=f[variance_i[(self.e1,self.e2)]].data
        
            self.wcs=pywcs.WCS(e.header)


            [e1,e2],self.fluxes,self.flux_errors=self.fit_image(rate,variance)

            spectra.append(concatenate(([self.e1,self.e2],self.fluxes,self.flux_errors)))
        
        savetxt("spectra.txt",spectra)
        self.spectra=spectra

