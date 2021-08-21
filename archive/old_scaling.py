# THIS FUNCTION IS NO LONGER USED
def power_scale(b_map, p, d_range):
	"""Takes as inputs a surface brightness map in Jy/str, a value for the power scaling and a data_range (list)
	   First, data is clipped so it only lies within the given data range. It is then scaled to lie between 0 and
	   and 10^p where p is the scale. Then take log of these values, giving a range of 1-p. These values will be 
	   scaled linearly to the set of availabel colors. (NOTE: This assumes negative values of p, can also take 
	   positive values, not described here"""
	
	data=b_map.value #[0,0,:,:] only keep dimensions we want


	if p < 0:
		power=-1.*p
		scaled_data=((data-np.min(data))/np.max(data))*(10**power)
		#print(np.min(scaled_data), np.max(scaled_data), scaled_data)
		output=np.log10(scaled_data+1)
		#print(np.min(output), np.max(output), output)
		return output
	elif p >0:
		scaled_data=((data-np.min(data))/np.max(data))*(p)
		output=10**scaled_data
		return output
	else:
		return data
