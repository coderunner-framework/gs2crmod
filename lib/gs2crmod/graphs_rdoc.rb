class CodeRunner::Gs2::GraphKits
# Graph of apar^2 vs time integrated over all space. No options
def apar2_vs_time_graphkit
end
# Plot the eigenfunction along the extended domain. Options mag, norm, z can be specified by using a short hand in the name of the graph, eg. efnmagnormz, efnmag, efnnorm etc. If the range is set to 0, it plots the whole eigenfunction. Otherwise it plot a small bit of it. Only specify kx or kx_index if magnetic shear is 0.
# Options:
#
# mag: Plot the magnitude, e.g. mag: true
#
# norm: Normalise the graph so that its maximum is 1, e.g. norm: true
#
# z: Plot quantities vs z = theta/shat rather than theta. See Beer, Cowley Hammet 1996, eg. z: true
#
# flip: Flip the y axis,  e.g. flip: true
#
# range: 
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# strongest_non_zonal_mode: Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true
def efn_graphkit
end
# Plot the eigenfunction along the extended domain. Options mag, norm, z can be specified by using a short hand in the name of the graph, eg. efnmagnormz, efnmag, efnnorm etc. If the range is set to 0, it plots the whole eigenfunction. Otherwise it plot a small bit of it. Only specify kx or kx_index if magnetic shear is 0.
# Options:
#
# mag: Plot the magnitude, e.g. mag: true
#
# norm: Normalise the graph so that its maximum is 1, e.g. norm: true
#
# z: Plot quantities vs z = theta/shat rather than theta. See Beer, Cowley Hammet 1996, eg. z: true
#
# flip: Flip the y axis,  e.g. flip: true
#
# range: 
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# strongest_non_zonal_mode: Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true
def efnim_graphkit
end
# Plot the eigenfunction along the extended domain. Options mag, norm, z can be specified by using a short hand in the name of the graph, eg. efnmagnormz, efnmag, efnnorm etc. If the range is set to 0, it plots the whole eigenfunction. Otherwise it plot a small bit of it. Only specify kx or kx_index if magnetic shear is 0.
# Options:
#
# mag: Plot the magnitude, e.g. mag: true
#
# norm: Normalise the graph so that its maximum is 1, e.g. norm: true
#
# z: Plot quantities vs z = theta/shat rather than theta. See Beer, Cowley Hammet 1996, eg. z: true
#
# flip: Flip the y axis,  e.g. flip: true
#
# range: 
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# strongest_non_zonal_mode: Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true
def efnmag_graphkit
end
# Plot the eigenfunction along the extended domain. Options mag, norm, z can be specified by using a short hand in the name of the graph, eg. efnmagnormz, efnmag, efnnorm etc. If the range is set to 0, it plots the whole eigenfunction. Otherwise it plot a small bit of it. Only specify kx or kx_index if magnetic shear is 0.
# Options:
#
# mag: Plot the magnitude, e.g. mag: true
#
# norm: Normalise the graph so that its maximum is 1, e.g. norm: true
#
# z: Plot quantities vs z = theta/shat rather than theta. See Beer, Cowley Hammet 1996, eg. z: true
#
# flip: Flip the y axis,  e.g. flip: true
#
# range: 
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# strongest_non_zonal_mode: Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true
def eigenfunction_graphkit
end
# 'es_heat_by_ky_vs_time' or 'es_heat_by_kx_vs_time': Electrostatic Heat Flux vs Time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_heat_by_kx_vs_time_graphkit
end
# 'es_heat_by_ky_vs_time' or 'es_heat_by_kx_vs_time': Electrostatic Heat Flux vs Time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_heat_by_kxy_or_mode_vs_time_graphkit
end
# 'es_heat_by_ky_vs_time' or 'es_heat_by_kx_vs_time': Electrostatic Heat Flux vs Time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_heat_by_ky_vs_time_graphkit
end
# 'es_heat_by_ky_vs_time' or 'es_heat_by_kx_vs_time': Electrostatic Heat Flux vs Time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_heat_by_mode_vs_time_graphkit
end
# Graph of electrostatic contribution to heat flux at a given time vs kx and ky
# Options:
#
# with: Gnuplot Option (may not apply when using other packages), e.g. with: 'lp' or with 'pm3d palette'
def es_heat_flux_vs_ky_vs_kx_graphkit
end
# Heat flux vs time for each species.
# Options:
#
# t_index_window: [begin, end], window of time indices to plot (e.g. t_index_window: [0,10])
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_heat_flux_vs_time_graphkit
end
# Momentum flux vs time for each species.
# Options:
#
# t_index_window: [begin, end], window of time indices to plot (e.g. t_index_window: [0,10])
#
# species_index: Which GS2 species to plot the graph for (1-based).
def es_mom_flux_vs_time_graphkit
end
# 'growth_rate_by_ky_vs_time' or 'growth_rate_by_kx_vs_time': Growth rate vs time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def growth_rate_by_kx_vs_time_graphkit
end
# 'growth_rate_by_ky_vs_time' or 'growth_rate_by_kx_vs_time': Growth rate vs time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def growth_rate_by_kxy_or_mode_vs_time_graphkit
end
# 'growth_rate_by_ky_vs_time' or 'growth_rate_by_kx_vs_time': Growth rate vs time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def growth_rate_by_ky_vs_time_graphkit
end
# 'growth_rate_by_ky_vs_time' or 'growth_rate_by_kx_vs_time': Growth rate vs time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def growth_rate_by_mode_vs_time_graphkit
end
# growth_rate_vs_ky or growth_rate_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /growth_rate_vs_kx_vs_ky/. 
def growth_rate_vs_kx_graphkit
end
# 3D plot of growth rates vs ky and kx for phi^2
# Options:
#
# with: Gnuplot Option (may not apply when using other packages), e.g. with: 'lp' or with 'pm3d palette'
def growth_rate_vs_kx_vs_ky_graphkit
end
# growth_rate_vs_ky or growth_rate_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /growth_rate_vs_kx_vs_ky/. 
def growth_rate_vs_kxy_graphkit
end
# growth_rate_vs_ky or growth_rate_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /growth_rate_vs_kx_vs_ky/. 
def growth_rate_vs_ky_graphkit
end
# Graph of total heat flux vs time. No options
def hflux_tot_vs_time_graphkit
end
# Graph of the k_parallel at a given kx and ky
# Options:
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# strongest_non_zonal_mode: Plot the graph requested for the mode with the highest value of phi^2. Overrides ky, kx, ky_index, kx_index. Can be set true or false; e.g. strongest_non_zonal_mode: true
def kpar_spectrum_graphkit
end
# ky_spectrum or kx_spectrum: Graph of phi^2 vs kx or ky
# Options:
#
# t: float, value of time at which to plot (e.g. t: 2.45)
#
# t_index: The (1-based) time index
def kx_spectrum_graphkit
end
# ky_spectrum or kx_spectrum: Graph of phi^2 vs kx or ky
# Options:
#
# t: float, value of time at which to plot (e.g. t: 2.45)
#
# t_index: The (1-based) time index
def kxy_spectrum_graphkit
end
# ky_spectrum or kx_spectrum: Graph of phi^2 vs kx or ky
# Options:
#
# t: float, value of time at which to plot (e.g. t: 2.45)
#
# t_index: The (1-based) time index
def ky_spectrum_graphkit
end
# A graph of the evolution of a single Lagrangian kx vs Eulerian kx and ky. Principally for debugging purposes
def lagrangian_kx_graphkit
end
# The potential at the outboard midplane
# Options:
#
# rgbformulae: Gnuplot Option (may not apply when using other packages), sets colour mapping. See gnuplot help set rgbformulae
#
# limit: Limit the range of quantity begin plotted - any values of the quantity outside the limits will be set to the limit: eg. limit: [0,80]
def phi0_vs_x_vs_y_graphkit
end
# 'phi2_by_ky_vs_time' or 'phi2_by_kx_vs_time': Phi^2 over time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def phi2_by_kx_vs_time_graphkit
end
# 'phi2_by_ky_vs_time' or 'phi2_by_kx_vs_time': Phi^2 over time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def phi2_by_kxy_or_mode_vs_time_graphkit
end
# 'phi2_by_ky_vs_time' or 'phi2_by_kx_vs_time': Phi^2 over time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def phi2_by_ky_vs_time_graphkit
end
# 'phi2_by_ky_vs_time' or 'phi2_by_kx_vs_time': Phi^2 over time for a given kx or ky, integrated over the other direction
# Options:
#
# ky: float, value of ky at which to plot (e.g. ky: 0.1)
#
# ky_index: integer, index of ky at which to plot (e.g. ky_index: 20)
#
# kx: float, value of kx at which to plot (e.g. kx: 0.1)
#
# kx_index: integer, index of kx at which to plot (e.g. kx_index: 20)
def phi2_by_mode_vs_time_graphkit
end
# Graph of phi^2 vs time integrated over all space. No options
def phi2tot_vs_time_graphkit
end
# The potential as a function of cartesian coordinates, on one specified side of the flux tube (specified using the options :coordinate (:x, :y, :theta) and :side (:min, :max))
# Options:
#
# Rgeo: 
#
# n0:  The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is
#
# rho_star:  The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity
#
# t_index: The (1-based) time index
#
# nakx: The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4
#
# naky: The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4
#
# gs2_coordinate_factor: When set to 1, plot the graph in GS2 coordinates. When set to  0 plot the graph in real space. Can be set at any value between 0 and 1: the graph will smoothly distort between the two limits
#
# xmax: The (0-based) index of the maximum value of x to include in the plot
#
# xmin: The (0-based) index of the minimum value of x to include in the plot
#
# ymax: The (0-based) index of the maximum value of y to include in the plot
#
# ymin: The (0-based) index of the minimum value of y to include in the plot
#
# thetamax: The (0-based) index of the maximum value of theta to include in the plot
#
# thetamin: The (0-based) index of the minimum value of theta to include in the plot
#
# ncopies:  The number of periodic copies of the flux tube to include
#
# side: 
#
# coordinate: 
#
# torphi_values: An array of two values of the toroidal angle. The graph will be plotted in between those two values with poloidal cross sections at either end
def phi_flux_tube_boundary_surface_graphkit
end
# The potential as a function of y, x and theta
# Options:
#
# rgbformulae: Gnuplot Option (may not apply when using other packages), sets colour mapping. See gnuplot help set rgbformulae
#
# limit: Limit the range of quantity begin plotted - any values of the quantity outside the limits will be set to the limit: eg. limit: [0,80]
#
# t_index: The (1-based) time index
def phi_gs2_space_graphkit
end
# The potential as a function of cartesian coordinates, plus a close up section.
# Options:
#
# Rgeo: 
#
# n0:  The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is
#
# rho_star:  The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity
#
# t_index: The (1-based) time index
#
# nakx: The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4
#
# naky: The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4
#
# xmax: The (0-based) index of the maximum value of x to include in the plot
#
# xmin: The (0-based) index of the minimum value of x to include in the plot
#
# thetamax: The (0-based) index of the maximum value of theta to include in the plot
#
# thetamin: The (0-based) index of the minimum value of theta to include in the plot
#
# magnify:  The magnification factor of the small section. It can take any value greater than or equal to 1
def phi_magnifying_glass_graphkit
end
# The potential as a function of cartesian coordinates
# Options:
#
# rgbformulae: Gnuplot Option (may not apply when using other packages), sets colour mapping. See gnuplot help set rgbformulae
#
# limit: Limit the range of quantity begin plotted - any values of the quantity outside the limits will be set to the limit: eg. limit: [0,80]
#
# t_index: The (1-based) time index
def phi_real_space_graphkit
end
# The potential as a function of cartesian coordinates showing a cut at one toroidal angle, with multiple periodic copies of the flux tube used to fill the whole circle..
# Options:
#
# Rgeo: 
#
# n0:  The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is
#
# rho_star:  The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity
#
# t_index: The (1-based) time index
#
# nakx: The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4
#
# naky: The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4
#
# xmax: The (0-based) index of the maximum value of x to include in the plot
#
# xmin: The (0-based) index of the minimum value of x to include in the plot
#
# thetamax: The (0-based) index of the maximum value of theta to include in the plot
#
# thetamin: The (0-based) index of the minimum value of theta to include in the plot
#
# torphi: 
def phi_real_space_poloidal_plane_graphkit
end
# The potential as a function of cartesian coordinates showing showing the standard way of representing the turbulence, with two poloidal cuts and the inner and outer radial surfaces. Multiple copies of the flux tube are used to fill the space.
# Options:
#
# Rgeo: 
#
# n0:  The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is
#
# rho_star:  The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity
#
# t_index: The (1-based) time index
#
# nakx: The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4
#
# naky: The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4
#
# xmax: The (0-based) index of the maximum value of x to include in the plot
#
# xmin: The (0-based) index of the minimum value of x to include in the plot
#
# thetamax: The (0-based) index of the maximum value of theta to include in the plot
#
# thetamin: The (0-based) index of the minimum value of theta to include in the plot
#
# torphi_values: An array of two values of the toroidal angle. The graph will be plotted in between those two values with poloidal cross sections at either end
def phi_real_space_standard_representation_graphkit
end
# The potential as a function of cartesian coordinates, plotted on the six outer surfaces of constant x, y and theta.
# Options:
#
# Rgeo: 
#
# n0:  The toroidal mode number of the longest y mode. In effect it is the number of periodic copies of the flux tube that will fit in the torus. Periodicity requires that n0 q  is also an integer. If you specify :n0 where this is not the case, q will automatically be adjusted until it is
#
# rho_star:  The ratio of the reference Lamour radius to the GS2 normalising length a. Cannot be specified at the same time as n0. If specified, both n0 and q will be adjusted to ensure periodicity
#
# t_index: The (1-based) time index
#
# nakx: The number of radial wave numbers to include in the plot. In effect, it is a low pass filter which reduces the resolution in the radial direction without changing the shape of the final surface. Minimum value is 4
#
# naky: The number of kys to include in the plot. In effect, it is a low pass filter which reduces the resolution in the y direction without changing the shape of the final surface. Minimum value is 4
#
# gs2_coordinate_factor: When set to 1, plot the graph in GS2 coordinates. When set to  0 plot the graph in real space. Can be set at any value between 0 and 1: the graph will smoothly distort between the two limits
#
# xmax: The (0-based) index of the maximum value of x to include in the plot
#
# xmin: The (0-based) index of the minimum value of x to include in the plot
#
# ymax: The (0-based) index of the maximum value of y to include in the plot
#
# ymin: The (0-based) index of the minimum value of y to include in the plot
#
# thetamax: The (0-based) index of the maximum value of theta to include in the plot
#
# thetamin: The (0-based) index of the minimum value of theta to include in the plot
#
# ncopies:  The number of periodic copies of the flux tube to include
def phi_real_space_surface_graphkit
end
# Graph of phi^2 at a given time vs kx and ky
# Options:
#
# with: Gnuplot Option (may not apply when using other packages), e.g. with: 'lp' or with 'pm3d palette'
def spectrum_graphkit
end
# Graph of phi^2 * ky^2 at a given time vs kpar and ky
# Options:
#
# with: Gnuplot Option (may not apply when using other packages), e.g. with: 'lp' or with 'pm3d palette'
#
# log: Plot the log of a given quantity (exact meaning varies). boolean
#
# no_zonal: Don't plot the ky=0 part (boolean, e.g. no_zonal: true)
#
# no_kpar0: Don't plot the kpar=0 part (boolean, e.g. no_kpar0: true)
def spectrum_vs_kpar_vs_ky_graphkit
end
# transient_amplification_vs_ky or transient_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_amplification_vs_kx_vs_ky/. 
def transient_amplification_vs_kx_graphkit
end
# transient_amplification_vs_ky or transient_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_amplification_vs_kx_vs_ky/. 
def transient_amplification_vs_kxy_graphkit
end
# transient_amplification_vs_ky or transient_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_amplification_vs_kx_vs_ky/. 
def transient_amplification_vs_ky_graphkit
end
# transient_es_heat_flux_amplification_vs_ky or transient_es_heat_flux_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_es_heat_flux_amplification_vs_kx_vs_ky/. 
def transient_es_heat_flux_amplification_vs_kx_graphkit
end
# transient_es_heat_flux_amplification_vs_ky or transient_es_heat_flux_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_es_heat_flux_amplification_vs_kx_vs_ky/. 
def transient_es_heat_flux_amplification_vs_kxy_graphkit
end
# transient_es_heat_flux_amplification_vs_ky or transient_es_heat_flux_amplification_vs_kx. Growth rates vs either ky or kx for phi^2 integrated over the other direction. For growth rates at a specific kx AND ky, see /transient_es_heat_flux_amplification_vs_kx_vs_ky/. 
def transient_es_heat_flux_amplification_vs_ky_graphkit
end
# Plots vspace diagnostics. All lines shouldn't stray much above 0.1 - otherwise large amounts of the distribution function is in the higher k velocity space and velocity space is probably unresolved. (NB This graph is here temporarily (ha ha) until I add the vspace diagnostics to the NetCDF file (or the apocalypse, whichever is sooner) EGH)
def vspace_diagnostics_graphkit
end
# zonal_spectrum: Graph of kx^4 phi^2 vs kx for ky=0
# Options:
#
# t: float, value of time at which to plot (e.g. t: 2.45)
#
# t_index: The (1-based) time index
def zonal_spectrum_graphkit
end
end
