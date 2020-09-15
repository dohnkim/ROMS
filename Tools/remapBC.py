import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date

import pyroms
import pyroms_toolbox

class nctime(object):
    pass

def remapBC(src_file, src_varname, wts_file, src_grd, dst_grd, dst_file, dxy=20, cdepth=0, kk=2, verbose=True):

    if src_varname == 'surf_el':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        dst_varname_north = 'zeta_north'
        dimensions_north = ('ocean_time', 'xi_rho')
        long_name_north = 'free-surface north boundary condition'
        field_north = 'zeta_north, scalar, series'
        dst_varname_south = 'zeta_south'
        dimensions_south = ('ocean_time', 'xi_rho')
        long_name_south = 'free-surface south boundary condition'
        field_south = 'zeta_south, scalar, series'
        dst_varname_east = 'zeta_east'
        dimensions_east = ('ocean_time', 'eta_rho')
        long_name_east = 'free-surface east boundary condition'
        field_east = 'zeta_east, scalar, series'
        dst_varname_west = 'zeta_west'
        dimensions_west = ('ocean_time', 'eta_rho')
        long_name_west = 'free-surface west boundary condition'
        field_west = 'zeta_west, scalar, series'
        units = 'meter'
    elif src_varname == 'water_temp':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'temp'
        dst_varname_north = 'temp_north'
        dimensions_north = ('ocean_time', 's_rho', 'xi_rho')
        long_name_north = 'potential temperature north boundary condition'
        field_north = 'temp_north, scalar, series'
        dst_varname_south = 'temp_south'
        dimensions_south = ('ocean_time', 's_rho', 'xi_rho')
        long_name_south = 'potential temperature south boundary condition'
        field_south = 'temp_south, scalar, series'
        dst_varname_east = 'temp_east'
        dimensions_east = ('ocean_time', 's_rho', 'eta_rho')
        long_name_east = 'potential temperature east boundary condition'
        field_east = 'temp_east, scalar, series'
        dst_varname_west = 'temp_west'
        dimensions_west = ('ocean_time', 's_rho', 'eta_rho')
        long_name_west = 'potential temperature west boundary condition'
        field_west = 'temp_west, scalar, series'
        units = 'Celsius'        
    elif src_varname == 'salinity':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'salt'
        dst_varname_north = 'salt_north'
        dimensions_north = ('ocean_time', 's_rho', 'xi_rho')
        long_name_north = 'salinity north boundary condition'
        field_north = 'salt_north, scalar, series'
        dst_varname_south = 'salt_south'
        dimensions_south = ('ocean_time', 's_rho', 'xi_rho')
        long_name_south = 'salinity south boundary condition'
        field_south = 'salt_south, scalar, series'
        dst_varname_east = 'salt_east'
        dimensions_east = ('ocean_time', 's_rho', 'eta_rho')
        long_name_east = 'salinity east boundary condition'
        field_east = 'salt_east, scalar, series'
        dst_varname_west = 'salt_west'
        dimensions_west = ('ocean_time', 's_rho', 'eta_rho')
        long_name_west = 'salinity west boundary condition'
        field_west = 'salt_west, scalar, series'
        units = 'PSU'        
    else:
        raise ValueError('Undefined src_varname')
    
    src_nc = Dataset(src_file) 
    # read time and change units to 'days' for ROMS
    dtime = num2date(src_nc.variables['time'][0],units=src_nc['time'].units,calendar=src_nc['time'].calendar)
    time  = date2num(dtime,units='days since 2000-01-01 00:00:00',calendar='gregorian')
    
    # get time
    nctime.long_name = 'time'
    nctime.units = 'days' # for ROMS
    #nctime.calendar = 'gregorian'

    # create destination file
    print('==>Info: Creating file', dst_file)
    if os.path.exists(dst_file) is True: os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_bdry_file(dst_file, dst_grd, nctime)

    # open BC file
    dst_nc = Dataset(dst_file, 'a', format='NETCDF3_64BIT')

    #load var & get missing value
    spval   = src_nc[src_varname]._FillValue
    src_var = src_nc[src_varname][0] # remove time

    # determine variable dimension
    #ndim = len(src_var.dimensions)
    ndim = len(src_var.shape)

    if ndim == 3:
        # build intermediate zgrid
        zlevel = -z[::-1,0,0]
        nzlevel = len(zlevel)
        dst_zcoord = pyroms.vgrid.z_coordinate(dst_grd.vgrid.h, zlevel, nzlevel)
        dst_grdz = pyroms.grid.ROMS_Grid(dst_file+'_Z', dst_grd.hgrid, dst_zcoord)

        src_var = src_var[:,1:-1,1:-1] # remove 4 corners
    else: src_var = src_var[1:-1,1:-1] # remove 4 corners
    
    # create variable in boudary file
    if verbose: print('==>Info: Creating variable', dst_varname_north)
    _=dst_nc.createVariable(dst_varname_north, 'f8', dimensions_north, fill_value=spval)
    dst_nc.variables[dst_varname_north].long_name = long_name_north
    dst_nc.variables[dst_varname_north].units = units
    dst_nc.variables[dst_varname_north].field = field_north

    if verbose: print('==>Info: Creating variable', dst_varname_south)
    _=dst_nc.createVariable(dst_varname_south, 'f8', dimensions_south, fill_value=spval)
    dst_nc.variables[dst_varname_south].long_name = long_name_south
    dst_nc.variables[dst_varname_south].units = units
    dst_nc.variables[dst_varname_south].field = field_south

    if verbose: print('==>Info: Creating variable', dst_varname_east)
    _=dst_nc.createVariable(dst_varname_east, 'f8', dimensions_east, fill_value=spval)
    dst_nc.variables[dst_varname_east].long_name = long_name_east
    dst_nc.variables[dst_varname_east].units = units
    dst_nc.variables[dst_varname_east].field = field_east

    if verbose: print('==>Info: Creating variable', dst_varname_west)
    _=dst_nc.createVariable(dst_varname_west, 'f8', dimensions_west, fill_value=spval)
    dst_nc.variables[dst_varname_west].long_name = long_name_west
    dst_nc.variables[dst_varname_west].units = units
    dst_nc.variables[dst_varname_west].field = field_west
    
    # remapping
    if verbose: print('==>Info: remapping', dst_varname, ', time =', time)


    if ndim == 3:
        # flood the grid
        if verbose: print('==>Info: flood the grid.', src_var.shape)
        src_varz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_var, src_grd, pos=pos, spval=spval, \
                                dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    if verbose: print('==>Info: horizontal interpolation using scrip weights\n',
                      '         about to call remap %s\n'%wts_file,
                      '         ', src_varz.shape)
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        if verbose: print('==>Info: vertical interpolation from standard z level to sigma')
        dst_var_north = pyroms.remapping.z2roms(dst_varz[::-1, Mp-1:Mp, :], \
                                                dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                                                flood=False, irange=(0,Lp), jrange=(Mp-1,Mp))
        dst_var_south = pyroms.remapping.z2roms(dst_varz[::-1, 0:1, :], \
                                                dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                                                flood=False, irange=(0,Lp), jrange=(0,1))
        dst_var_east  = pyroms.remapping.z2roms(dst_varz[::-1, :, Lp-1:Lp], \
                                                dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                                                flood=False, irange=(Lp-1,Lp), jrange=(0,Mp))
        dst_var_west  = pyroms.remapping.z2roms(dst_varz[::-1, :, 0:1], \
                                                dst_grdz, dst_grd, Cpos=Cpos, spval=spval, \
                                                flood=False, irange=(0,1), jrange=(0,Mp))
    else:
        dst_var_north = dst_varz[-1, :]
        dst_var_south = dst_varz[ 0, :]
        dst_var_east  = dst_varz[ :,-1]
        dst_var_west  = dst_varz[ :, 0]

    # write data in destination file
    if verbose: print('==>Info: write data in destination file')
    dst_nc.variables['ocean_time'][0] = time
    dst_nc.variables[dst_varname_north][0] = np.squeeze(dst_var_north)
    dst_nc.variables[dst_varname_south][0] = np.squeeze(dst_var_south)
    dst_nc.variables[dst_varname_east][0]  = np.squeeze(dst_var_east)
    dst_nc.variables[dst_varname_west][0]  = np.squeeze(dst_var_west)

    # close destination file
    src_nc.close()
    dst_nc.close()

    if src_varname == 'surf_el':
        return dst_varz
    
