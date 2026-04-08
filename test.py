# 1 pixel = 0.267112299e-3 m in the experiments
# 1pixel = 1/50 s

def compute_u_vortex_exp(pixel_vertical, pixel_temps):
	"""Compute U_vortex from experimental data."""
	# Implementation for experimental computation
	return pixel_vertical * 0.267112299e-3 / (pixel_temps / 50)

print(compute_u_vortex_exp(101, 9))