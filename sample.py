import numpy as np
import MEX as MEX 
import scipy.interpolate as interpolate

loadhalopath = "/Users/Scott/Desktop/cosmo/halo_with_accelerations_no_disk.dat"
mass, r, theta_pos, phi_pos = MEX.LoadHaloalt(loadhalopath)
print(min(r))
l = [0, 2, 4, 6, 8, 10]
binss = [200, 250, 300, 350, 400]
for l_max in l: 
	print("l_max = ",l_max)
	for bins in binss : 
		print("bins =",bins)
		bincenters, r_trunc, log_bins = MEX.bin_r(r,bins,210)
		
		basis_functions = MEX.basis_function_mod(l_max, bincenters, log_bins,r, r_trunc, mass, theta_pos, phi_pos)
		
		int_lm1, int_lm2 = MEX.MEX_integral_lm_force(l_max,basis_functions,bincenters)
		
		# forces over array bincenter[1:]
		F_x = []
		F_y = []
		F_z = []

		i = 0 
		while i < len(theta_pos):

			f_x,f_y,f_z,f_r,f_phi,f_theta = MEX.MEX_force(l_max, theta_pos[i], phi_pos[i], bincenters , int_lm1, int_lm2)
			
			
			f_x = interpolate.interp1d(bincenters[1:],f_x,  kind='linear', fill_value='extrapolate' )
			f_y = interpolate.interp1d(bincenters[1:],f_y,  kind='linear', fill_value='extrapolate' )
			f_z = interpolate.interp1d(bincenters[1:],f_z,  kind='linear', fill_value='extrapolate' )
			
			f_xx = f_x(r[i])
			f_yy = f_y(r[i])
			f_zz = f_z(r[i])
			

			F_x.append(f_xx)
			F_y.append(f_yy)
			F_z.append(f_zz)
			

			i += 1




		name1 = '/Users/Scott/Desktop/cosmo/force/fx_%s_%s.out' %(l_max, bins)
		np.savetxt(name1, F_x, delimiter=',')

		name2 = '/Users/Scott/Desktop/cosmo/force/fy_%s_%s.out' %(l_max, bins)
		np.savetxt(name2, F_y, delimiter=',')

		name3 = '/Users/Scott/Desktop/cosmo/force/fz_%s_%s.out' %(l_max, bins)
		np.savetxt(name3, F_z, delimiter=',')









