import os
import numpy as np
from netCDF4 import Dataset, date2num, num2date

import pyroms
import pyroms_toolbox

class nctime(object):
    pass

def remapIC(src_file, src_varname, wts_file, src_grd, dst_grd, dst_file, dxy=20, cdepth=0, kk=0, verbose=True):

    if src_varname == 'surf_el':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'zeta'
        dimensions = ('ocean_time', 'eta_rho', 'xi_rho')
        long_name = 'free-surface'
        units = 'meter'
        field = 'free-surface, scalar, series'
    elif src_varname == 'water_temp':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'temp'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'potential temperature'
        units = 'Celsius'
        field = 'temperature, scalar, series'
    elif src_varname == 'salinity':
        pos = 't'; Cpos = 'rho'
        Mp, Lp = dst_grd.hgrid.mask_rho.shape; z = src_grd.z_t
        #wts_file = 'remap_weights_%s_to_%s_bilinear_t_to_rho.nc'%(src_name,dst_name)
        dst_varname = 'salt'
        dimensions = ('ocean_time', 's_rho', 'eta_rho', 'xi_rho')
        long_name = 'salinity'
        units = 'PSU'
        field = 'salinity, scalar, series'
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

    # create IC file
    print('\n==>Info: Creating file', dst_file)
    if os.path.exists(dst_file) is True: os.remove(dst_file)
    pyroms_toolbox.nc_create_roms_file(dst_file, dst_grd, nctime)

    # open IC file
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

    # create variable in file
    print('==>Info: Creating variable', dst_varname)
    _=dst_nc.createVariable(dst_varname, 'f8', dimensions, fill_value=spval)
    dst_nc.variables[dst_varname].long_name = long_name
    dst_nc.variables[dst_varname].units = units
    dst_nc.variables[dst_varname].field = field


    # remapping
    print('==>Info: remapping', dst_varname, ', time =', time)


    if ndim == 3:
        # flood the grid
        print('==>Info: flood the grid.', src_var.shape)
        src_varz = pyroms_toolbox.Grid_HYCOM.flood_fast(src_var, src_grd, pos=pos, spval=spval, \
                                dxy=dxy, cdepth=cdepth, kk=kk, verbose=verbose)
    else:
        src_varz = src_var

    # horizontal interpolation using scrip weights
    print('==>Info: horizontal interpolation using scrip weights')
    print('         about to call remap ' + wts_file)
    print('         ', src_varz.shape)
    dst_varz = pyroms.remapping.remap(src_varz, wts_file, spval=spval)

    if ndim == 3:
        # vertical interpolation from standard z level to sigma
        print('==>Info: vertical interpolation from standard z level to sigma')
        dst_var = pyroms.remapping.z2roms(dst_varz[::-1,:,:], dst_grdz, \
                          dst_grd, Cpos=Cpos, spval=spval, flood=False)
    else:
        dst_var = dst_varz

    # write data in destination file
    print('==>Info: write data in destination file')
    dst_nc.variables['ocean_time'][0] = time
    dst_nc.variables[dst_varname][0] = dst_var

    # close destination file
    src_nc.close()
    dst_nc.close()

    if src_varname == 'surf_el':
        return dst_varz
    
