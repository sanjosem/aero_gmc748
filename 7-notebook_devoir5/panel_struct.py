import numpy as np

class Panel:
    """
    Contains panel data structure
    """
    def __init__(self, x1, y1, x2, y2):
        """
        Sets the panel bound points and compute the panel center and length.
        
        Parameters
        ----------
        x1, y1 : coordinate of the first end point
        x2, y2 : coordinate of the second end point
        """
        from math import sqrt
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.xc = 0.5*(x1+x2)
        self.yc = 0.5*(y1+y2)
        self.length = sqrt((x2-x1)**2+(y2-y1)**2)

        self.vortex_strength = 0.
        self.source_strength = 0.

        # Angle between the panel and the (Ox) axis
        self.alpha = np.arctan2((y2-y1),(x2-x1))
        # Angle between the panel normal and the (Ox) axis
        self.beta = self.alpha + np.pi / 2.

    def source_influence_analytic(self,x,y,beta):
        """
        Compute the influence term at (x,y) for the source panel
        using the analytic expression 

        Parameters
        ----------
        x,y : numpy array
            coordinate of the colocation points where the influence is computed
        beta : numpy array (of the same size than x and y)
            angle measured from (Ox) of the direction on which the influence coefficient has to be projected

        Returns
        ----------
        I : influence parameter computed by integration
        """

        # global to panel coordinates : rotation by alpha of the reference frame
        xp = np.cos(self.alpha) * (x-self.x1) + np.sin(self.alpha) * (y-self.y1) 
        yp = - np.sin(self.alpha) * (x-self.x1) + np.cos(self.alpha) * (y-self.y1) 
        local_x1 = 0.
        local_x2 = self.length
        local_y1 = 0.
        local_y2 = 0.
        r1 = np.sqrt( xp**2 + yp**2 ) # local_x1=local_y1=0
        r2 = np.sqrt( (xp-local_x2)**2 + yp**2 ) # local_z2=0
        theta1 = np.arctan2( yp, xp )
        theta2 = np.arctan2( yp, (xp-local_x2) )
        up = np.log(r1/r2) / (2. * np.pi)
        vp = (theta2 - theta1) / (2. * np.pi)
        # panel velocity to global velocity : rotation by alpha of the vector components
        u =  np.cos(self.alpha) * up - np.sin(self.alpha) * vp
        v = np.sin(self.alpha) * up + np.cos(self.alpha) * vp
        # Influence = projection along provided beta angle (not self.beta)
        I = u * np.cos(beta) + v * np.sin(beta)
        return I

    def compute_source_analytic(self,X,Y):
        """
        Create the source panel and compute U,V,PHI,PSI

        Parameters
        ----------
        X,Y : numpy array 
            coordinate of the field points where the fields are computed

        """
        from elementary_flows import Distributed_constant_Source
        src=Distributed_constant_Source(self.x1,self.y1,self.x2,self.y2,self.source_strength)
        src.compute_velocity_analytic(X,Y)
        self.source = src

    
    


class NACA_Airfoil:
    """
    Contain NACA airfoil parameters and associated functions to compute the form
    """
    def __init__(self, number, chord,sharp=True):
        """
        Sets the airfoil parameter.
        
        Parameters
        ----------
        number: 4 digit string 
            NACA number
        chord : float
            chord of the airfoil
        sharp: boolean
            sharp or blunt TE (default sharp=True)

        """
        self.number = number
        self.camber_max = int(number[0])/100.
        self.camber_pos = int(number[1])/10.
        self.thickness = int(number[2:4])/100.
        self.chord = chord
        self.coeff = np.array((0.2966,-0.1260,-0.3516,0.2843,-0.1015))
        self.expo = np.array((0.5,1,2,3,4))
        if sharp:
            self.coeff[-1]=-0.1033

    def half_thickness(self,x):
        """
        Compute the half thickness of symmetric NACA 
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord
        def yt(xc):
            y_t = np.zeros_like(xc)
            for ic in range(5):
                y_t +=  self.coeff[ic] * xc**(self.expo[ic])
            y_t *= self.thickness / 0.2
            return y_t

        self.shape = yt(xc) * self.chord


    def camber(self,x):
        """
        Compute the mean camber line of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord

        def yc(xc):
            y_c = np.zeros_like(xc)
            xc_up = xc < self.camber_pos
            xc_dn = xc >= self.camber_pos
            if np.any(xc_up):
                y_c[xc_up] =  self.camber_max / self.camber_pos**2 * (2.*self.camber_pos*xc[xc_up]-xc[xc_up]**2)
            y_c[xc_dn] =  self.camber_max / (1.-self.camber_pos)**2 * (1.-2.*self.camber_pos + 2.*self.camber_pos*xc[xc_dn]-xc[xc_dn]**2)
            
            return y_c

        self.mean_camber_line = yc(xc) * self.chord

    def deriv_camber(self,x):
        """
        Compute the mean camber line derivative of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord
        """        
        xc = x / self.chord

        def dycdx(xc):
            dy_c = np.zeros_like(xc)
            xc_up = xc < self.camber_pos
            xc_dn = xc >= self.camber_pos
            if np.any(xc_up):
                dy_c[xc_up] =  2. * self.camber_max / self.camber_pos**2 * (self.camber_pos-xc[xc_up])
            dy_c[xc_dn] =  2.*self.camber_max / (1.-self.camber_pos)**2 * (self.camber_pos -xc[xc_dn])
            
            return dy_c

        self.deriv_camber_fct = dycdx(xc)

    def profile_side_shape(self,x):
        """
        Compute the suction side of NACA 4 digits
        
        Parameters
        ----------
        x: coordinates along the chord

        Returns
        -------
        xU: coordinates along the x-axis of the upper edge of the profile
        yU: coordinates along the y-axis of the upper edge of the profile
        xL: coordinates along the x-axis of the lower edge of the profile
        yL: coordinates along the y-axis of the lower edge of the profile
        """
        self.half_thickness(x)
        if self.camber_max>0:
            self.camber(x)
            self.deriv_camber(x)
            theta = np.arctan(self.deriv_camber_fct)
            xU = x - self.shape * np.sin(theta)
            yU = self.mean_camber_line + self.shape * np.cos(theta)
            xL = x + self.shape * np.sin(theta)
            yL = self.mean_camber_line - self.shape * np.cos(theta)
        else:
            xU = x
            yU = self.shape
            xL = x
            yL = - self.shape
        
        return xU,yU,xL,yL

