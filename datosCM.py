import numpy as  np
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
import gala.coordinates as gc

def datos_CM(cte):
	#datos de la LMC Sacados de van der Marel, ApJL 832, L23 (2016)
	#Tabla 2 columna (2) PMSTGAS
	ra_LMC=78.76*np.pi/180
	dec_LMC=-69.19*np.pi/180
	D_LMC=50.1 #en kpc
	vl_LMC=262.2 #en km/s
	mual_LMC=1.872 #en mas/yr
	mude_LMC=0.224#en mas/yr
	#Datosde SMC
	ra_SMC=12.80*np.pi/180
	dec_SMC=-73.15*np.pi/180
	D_SMC=62.80 #en kpc
	vl_SMC=145.6   #en km/s
	mual_SMC=0.79 #en mas/yr
	mude_SMC=-1.256 #en mas/yr
	#Seteamos los valores para el sistema galactico ya que los que vienen por default estan desactualizados
	v_sun = [-12.9, 245.6, 7.78] * (u.km / u.s)  # [vx, vy, vz]
	gc_frame = coord.Galactocentric(galcen_distance=8.122*u.kpc,galcen_v_sun=v_sun,z_sun=20.8*u.pc)
	#inicializamos el objeto skycoord y hacemos la corrección por el mov solar de la LMC
	LMC_car=SkyCoord(ra=ra_LMC*u.rad,dec=dec_LMC*u.rad,distance=D_LMC*u.kpc,frame='icrs')
	c_LMC=coord.SkyCoord(ra=ra_LMC*u.rad,dec=dec_LMC*u.rad,distance=D_LMC*u.kpc,pm_ra_cosdec=mual_LMC*u.mas/u.yr,pm_dec=mude_LMC*u.mas/u.yr,radial_velocity=vl_LMC*u.km/u.s,frame='icrs')
	corr_LMC=gc.reflex_correct(c_LMC,galactocentric_frame = gc_frame)
	mualpha_LMC=corr_LMC.pm_ra_cosdec.value
	mudelta_LMC=corr_LMC.pm_dec.value
	vlos_LMC=corr_LMC.radial_velocity.value
	#inicializamos el objeto skycoord y hacemos la corrección por el mov solar de la SMC
	SMC_car=SkyCoord(ra=ra_SMC*u.rad,dec=dec_SMC*u.rad,distance=D_SMC*u.kpc,frame='icrs')
	c_SMC=coord.SkyCoord(ra=ra_SMC*u.rad,dec=dec_SMC*u.rad,distance=D_SMC*u.kpc,pm_ra_cosdec=mual_SMC*u.mas/u.yr,pm_dec=mude_SMC*u.mas/u.yr,radial_velocity=vl_SMC*u.km/u.s,frame='icrs')
	corr_SMC=gc.reflex_correct(c_SMC,galactocentric_frame = gc_frame)
	mualpha_SMC=corr_SMC.pm_ra_cosdec.value
	mudelta_SMC=corr_SMC.pm_dec.value
	vlos_SMC=corr_SMC.radial_velocity.value
	#Velocidad del centro de masa
	#cte=9
	vcmx=(cte*corr_LMC.velocity.d_x.value+corr_SMC.velocity.d_x.value)/(cte+1)
	vcmy=(cte*corr_LMC.velocity.d_y.value+corr_SMC.velocity.d_y.value)/(cte+1)
	vcmz=(cte*corr_LMC.velocity.d_z.value+corr_SMC.velocity.d_z.value)/(cte+1)
	vsh=np.sqrt(vcmx**2+vcmy**2+vcmz**2)
	print('velocidad del CM',vcmx*u.km/u.s,vcmy*u.km/u.s,vcmz*u.km/u.s)
	print('modulo de velocidad',vsh*u.km/u.s)
	# posicion del centro de masa
	rcm=(cte*LMC_car.cartesian+SMC_car.cartesian)/(cte+1)
	return vsh,rcm,vcmx,vcmy,vcmz